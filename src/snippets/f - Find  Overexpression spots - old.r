

if (F)
{


        ##### Overexpression Cluster ######
  #####      old version      ######






  if (! ("summary.spot.cutoff" %in% names(preferences)))  preferences$summary.spot.cutoff = 0.90





  peaks = apply(apply(metadata, 2, function(x){ (x - min(x)) / (max(x) - min(x)) }), 1, max, na.rm=T)



  GS.infos.overexpression = list()

  GS.infos.overexpression$overview.map = peaks







  peaks = matrix(peaks, preferences$dim.som1, preferences$dim.som1)
  peaks[which(peaks < quantile(peaks, preferences$summary.spot.cutoff))] = NA

  e.cluster = matrix(NA, preferences$dim.som1, preferences$dim.som1)
  e.cluster[which(!is.na(peaks))] = -1
  count.cluster = 1


  while (any(!is.na(peaks)))
  {

    start.pix = which(peaks == max(peaks,na.rm=T), arr.ind=T)[1,]

    e.cluster = col.pix(e.cluster, start.pix[1], start.pix[2], count.cluster, preferences$dim.som1)

    peaks[which(e.cluster == count.cluster)] = NA

    count.cluster = count.cluster + 1
  }



  e.cluster = as.vector(e.cluster)

  for (i in 1:length(na.omit(unique(e.cluster))))
  {
    nodes = which(e.cluster == i)


    geneset.genes = rownames(indata)[which(som.nodes %in% nodes)]
    geneset.ids = unique(na.omit(gene.ids[geneset.genes]))


    GS.infos.overexpression[[LETTERS[i]]] = list()

    GS.infos.overexpression[[LETTERS[i]]]$metagenes = as.numeric(nodes)
    GS.infos.overexpression[[LETTERS[i]]]$genes = geneset.genes

    GS.infos.overexpression[[LETTERS[i]]]$mask = rep(NA, (preferences$dim.som1*preferences$dim.som1))
    GS.infos.overexpression[[LETTERS[i]]]$mask[as.numeric(nodes)] = 1

    if (preferences$geneset.analysis)
      GS.infos.overexpression[[LETTERS[i]]]$HG.p = GeneSet.HG(geneset.ids, unique.protein.ids, gs.def.list, sort=T)$p.value
  }



  o = order(sapply(GS.infos.overexpression[-1], function(x)
  {
    mean.spot.metagene = apply(metadata[x$metagenes, ,drop=F], 2, mean)
    return(which.max(mean.spot.metagene))
  }))
  o = c(1, o+1)

  GS.infos.overexpression = GS.infos.overexpression[o]
  names(GS.infos.overexpression)[-1] = LETTERS[1:(length(GS.infos.overexpression)-1)]



}





  ##### Overexpression Cluster ######
  ##### flooding hills version ######




  if (! ("summary.spot.minsize" %in% names(preferences)))  preferences$summary.spot.minsize = 4
  if (! ("summary.spot.maxsize" %in% names(preferences)))  preferences$summary.spot.maxsize = 10
  if (! ("summary.spot.startlevel" %in% names(preferences)))  preferences$summary.spot.startlevel = 0.9







  ##### Overexpression Spots ######

  peaks = apply(apply(metadata, 2, function(x){ (x - min(x)) / (max(x) - min(x)) }), 1, max)

  GS.infos.overexpression = list()
  GS.infos.overexpression$overview.map = peaks
  GS.infos.overexpression$overview.mask = rep(NA, preferences$dim.som1^2)
  GS.infos.overexpression$filtered = F
  GS.infos.overexpression$spots = list()


  peaks = matrix(peaks, preferences$dim.som1, preferences$dim.som1)

  count.spots = 1
  for (quant in c(seq(preferences$summary.spot.startlevel, 0.9, 0.01), seq(0.9,0.999,0.001), 0.9999, 0.99999))
  {

    peaks.copy = peaks
    peaks.copy[which(peaks.copy < quant)] = NA


    cluster.mask = matrix(NA, preferences$dim.som1, preferences$dim.som1)
    cluster.mask[which(!is.na(peaks.copy))] = -1
    count.cluster = 1
    while (any(!is.na(peaks.copy)))
    {
      start.pix = which(peaks.copy == max(peaks.copy,na.rm=T), arr.ind=T)[1,]
      cluster.mask = col.pix(cluster.mask, start.pix[1], start.pix[2], count.cluster, preferences$dim.som1)

      peaks.copy[which(cluster.mask == count.cluster)] = NA

      count.cluster = count.cluster + 1
    }


    if (quant == 0.99999)
    {
      nice.spots = na.omit(unique(as.vector(cluster.mask)))
    }
    else
    {
      nice.spots = which(table(cluster.mask) <= preferences$summary.spot.maxsize)
    }

    cluster.mask[which(!cluster.mask%in%nice.spots)] = NA
    peaks[which(cluster.mask%in%nice.spots)] = NA


    cluster.mask = as.vector(cluster.mask)
    for (i in nice.spots)
    {
      nodes = which(cluster.mask == i)

      geneset.genes = rownames(indata)[which(som.nodes %in% nodes)]

      if (length(nodes) >= preferences$summary.spot.minsize || quant == 0.99999)
      {
        GS.infos.overexpression$overview.mask[as.numeric(nodes)] = count.spots

        GS.infos.overexpression$spots[[LETTERS[count.spots]]] = list()

        GS.infos.overexpression$spots[[LETTERS[count.spots]]]$metagenes = as.numeric(nodes)
        GS.infos.overexpression$spots[[LETTERS[count.spots]]]$genes = geneset.genes

        GS.infos.overexpression$spots[[LETTERS[count.spots]]]$mask = rep(NA, (preferences$dim.som1*preferences$dim.som1))
        GS.infos.overexpression$spots[[LETTERS[count.spots]]]$mask[as.numeric(nodes)] = 1

        GS.infos.overexpression$spots[[LETTERS[count.spots]]]$position = apply(apply(som.result$code.sum[nodes, 1:2], 2, range), 2, mean) + 0.5

        GS.infos.overexpression$spots[[LETTERS[count.spots]]]$beta.score = get.beta.statistic(set.data=metadata[GS.infos.overexpression$spots[[LETTERS[count.spots]]]$metagenes,,drop=F])$beta.score

        count.spots = count.spots + 1
      }
    }

  }


  o = order(sapply(GS.infos.overexpression$spots, function(x)
  {
    mean.spot.metagene = apply(metadata[x$metagenes, ,drop=F], 2, mean)
    return(which.max(mean.spot.metagene))
  }))

  GS.infos.overexpression$spots = GS.infos.overexpression$spots[o]
  names(GS.infos.overexpression$spots) = LETTERS[1:length(GS.infos.overexpression$spots)]

  GS.infos.overexpression$overview.mask[!is.na(GS.infos.overexpression$overview.mask)] = match(GS.infos.overexpression$overview.mask[!is.na(GS.infos.overexpression$overview.mask)], o)


  GS.infos.overexpression$spotdata = t(sapply(GS.infos.overexpression$spots, function(x) if (length(x$genes > 0))  colMeans(indata[x$genes,,drop=F]) else rep(0, ncol(indata))))
  colnames(GS.infos.overexpression$spotdata) = colnames(indata)









  ##### Underexpression Spots ######


  peaks = apply(apply(metadata, 2, function(x){ (x - min(x)) / (max(x) - min(x)) }), 1, min)

  GS.infos.underexpression = list()
  GS.infos.underexpression$overview.map = peaks
  GS.infos.underexpression$overview.mask = rep(NA, preferences$dim.som1^2)
  GS.infos.underexpression$filtered = F
  GS.infos.underexpression$spots = list()


  peaks = matrix(peaks, preferences$dim.som1, preferences$dim.som1)

  count.spots = 1
  for (quant in c(seq(preferences$summary.spot.startlevel, 0.9, 0.01), seq(0.9,0.999,0.001), 0.9999, 0.99999))
  {

    peaks.copy = peaks
    peaks.copy[which(peaks.copy > 1-quant)] = NA


    cluster.mask = matrix(NA, preferences$dim.som1, preferences$dim.som1)
    cluster.mask[which(!is.na(peaks.copy))] = -1
    count.cluster = 1
    while (any(!is.na(peaks.copy)))
    {
      start.pix = which(peaks.copy == max(peaks.copy,na.rm=T), arr.ind=T)[1,]
      cluster.mask = col.pix(cluster.mask, start.pix[1], start.pix[2], count.cluster, preferences$dim.som1)

      peaks.copy[which(cluster.mask == count.cluster)] = NA

      count.cluster = count.cluster + 1
    }


    if (quant == 0.99999)
    {
      nice.spots = na.omit(unique(as.vector(cluster.mask)))
    }
    else
    {
      nice.spots = which(table(cluster.mask) <= preferences$summary.spot.maxsize)
    }

    cluster.mask[which(!cluster.mask%in%nice.spots)] = NA
    peaks[which(cluster.mask%in%nice.spots)] = NA


    cluster.mask = as.vector(cluster.mask)
    for (i in nice.spots)
    {
      nodes = which(cluster.mask == i)

      geneset.genes = rownames(indata)[which(som.nodes %in% nodes)]

      if (length(nodes) >= preferences$summary.spot.minsize || quant == 0.99999)
      {
        GS.infos.underexpression$overview.mask[as.numeric(nodes)] = count.spots

        GS.infos.underexpression$spots[[letters[count.spots]]] = list()

        GS.infos.underexpression$spots[[letters[count.spots]]]$metagenes = as.numeric(nodes)
        GS.infos.underexpression$spots[[letters[count.spots]]]$genes = geneset.genes

        GS.infos.underexpression$spots[[letters[count.spots]]]$mask = rep(NA, (preferences$dim.som1*preferences$dim.som1))
        GS.infos.underexpression$spots[[letters[count.spots]]]$mask[as.numeric(nodes)] = 1

        GS.infos.underexpression$spots[[letters[count.spots]]]$position = apply(apply(som.result$code.sum[nodes, 1:2], 2, range), 2, mean) + 0.5

        GS.infos.underexpression$spots[[letters[count.spots]]]$beta.score = get.beta.statistic(set.data=metadata[GS.infos.underexpression$spots[[letters[count.spots]]]$metagenes,,drop=F])$beta.score

        count.spots = count.spots + 1
      }
    }

  }


  o = order(sapply(GS.infos.underexpression$spots, function(x)
  {
    mean.spot.metagene = apply(metadata[x$metagenes, ,drop=F], 2, mean)
    return(which.min(mean.spot.metagene))
  }))

  GS.infos.underexpression$spots = GS.infos.underexpression$spots[o]
  names(GS.infos.underexpression$spots) = letters[1:length(GS.infos.underexpression$spots)]

  GS.infos.underexpression$overview.mask[!is.na(GS.infos.underexpression$overview.mask)] = match(GS.infos.underexpression$overview.mask[!is.na(GS.infos.underexpression$overview.mask)], o)


  GS.infos.underexpression$spotdata = t(sapply(GS.infos.underexpression$spots, function(x) if (length(x$genes > 0))  colMeans(indata[x$genes,,drop=F]) else rep(0, ncol(indata))))
  colnames(GS.infos.underexpression$spotdata) = colnames(indata)









