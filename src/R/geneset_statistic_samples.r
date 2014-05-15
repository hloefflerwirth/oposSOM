pipeline.genesetStatisticSamples <- function()
{
  progress.current <- 0
  progress.max <- ncol(indata)
  util.progress(progress.current, progress.max)

  ### init parallel computing ###
  cl <- makeCluster(preferences$max.parallel.cores)

  ### perform GS analysis ###
  t.ensID.m <<- t.g.m[which(rownames(indata) %in% names(gene.ids)),]
  t.ensID.m <<- do.call(rbind, by(t.ensID.m, gene.ids, colMeans))

  if (preferences$geneset.analysis.exact)
  {
    gs.null.list <- list()

    for (i in 1:length(gs.def.list))
    {
      gs.null.list[[i]] <-
        list(Genes=sample(unique.protein.ids, length(gs.def.list[[i]]$Genes)))
    }

    null.scores <- parSapply(cl, 1:ncol(indata), function(m)
    {
      all.gene.statistic <- t.ensID.m[,m]
      spot.gene.ids <- unique.protein.ids

      scores <- GeneSet.GSZ(spot.gene.ids, all.gene.statistic, gs.null.list)

      for (spot.i in 1:length(spot.list.samples[[m]]$spots))
      {
        spot.genes <- spot.list.samples[[m]]$spots[[spot.i]]$genes
        spot.gene.ids <- unique(na.omit(gene.ids[spot.genes]))
        all.gene.statistic <- t.ensID.m[, m]

        scores <-
          c(scores, GeneSet.GSZ(spot.gene.ids, all.gene.statistic, gs.null.list))
      }

      return(scores)
    })

    null.culdensity <- ecdf(abs(unlist(null.scores)))
  }

  progress.current <- progress.max / 2
  util.progress(progress.current, progress.max)

  spot.list.samples <<- parLapply(cl, 1:length(spot.list.samples) , function(m)
  {
    x <- spot.list.samples[[m]]
    all.gene.statistic <- t.ensID.m[,m]
    spot.gene.ids <- unique.protein.ids

    x$GSZ.score <-
      GeneSet.GSZ(spot.gene.ids, all.gene.statistic, gs.def.list, sort=F)

    if (preferences$geneset.analysis.exact)
    {
      x$GSZ.p.value <- 1 - null.culdensity(abs(x$GSZ.score))
      names(x$GSZ.p.value) <- names(x$GSZ.score)
    }

    for (spot.i in 1:length(x$spots))
    {
      spot.genes <- x$spots[[spot.i]]$genes
      spot.gene.ids <- unique(na.omit(gene.ids[spot.genes]))

      if (length(spot.gene.ids) > 0)
      {
        x$spots[[spot.i]]$GSZ.score <-
          GeneSet.GSZ(spot.gene.ids, all.gene.statistic, gs.def.list, sort=F)

        x$spots[[spot.i]]$Fisher.p <-
          GeneSet.Fisher(spot.gene.ids, unique.protein.ids, gs.def.list, sort=T)
      } else
      {
        x$spots[[spot.i]]$GSZ.score <- rep(0, length(gs.def.list))
        names(x$spots[[spot.i]]$GSZ.score) <- names(gs.def.list)
        x$spots[[spot.i]]$Fisher.p <- rep(1, length(gs.def.list))
        names(x$spots[[spot.i]]$Fisher.p) <- names(gs.def.list)
      }

      if (preferences$geneset.analysis.exact)
      {
        x$spots[[spot.i]]$GSZ.p.value <-
          1 - null.culdensity(abs(x$spots[[spot.i]]$GSZ.score))

        names(x$spots[[spot.i]]$GSZ.p.value) <- names(x$spots[[spot.i]]$GSZ.score)
      }
    }

    return(x)
  })

  progress.current <- progress.max * 0.9
  util.progress(progress.current, progress.max)

  ### stop parallel computing ###
  try({ stopCluster(cl) }, silent=T)
}
