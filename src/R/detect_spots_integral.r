pipeline.detectSpotsIntegral <- function()
{
  metadata.scaled <- apply(metadata, 2, function(x)
  {
    (x-min(x)) / (max(x)-min(x))
  })

  ##### Overexpression Spots ######

  # extract sample modules
  sample.spot.list <- list()
  sample.spot.core.list <- list()
  n.sample.modules <- 0

  for (m in 1:ncol(indata))
  {
    # define bigger core regions
    core <- matrix(NA, preferences$dim.1stLvlSom, preferences$dim.1stLvlSom)
    core[which(metadata.scaled[,m] > preferences$spot.threshold.modules)] <- -1

    spot.i = 0

    while (nrow(which(core == -1, arr.ind=T)) > 0)
    {
      start.pix <- which(core == -1, arr.ind=T)[1,]
      spot.i <- spot.i + 1
      core <- col.pix(core, start.pix[1], start.pix[2], spot.i, preferences$dim.1stLvlSom)
    }

    # shrink each separate region to core size
    for (s.i in 1:max(core, na.rm=T))
    {
      if (sum(core == s.i, na.rm=T) > preferences$spot.coresize.modules)
      {
        core.metagenes <- which(core == s.i)
        o <- order(metadata[core.metagenes,m], decreasing=T)[1:preferences$spot.coresize.modules]
        core[setdiff(core.metagenes, core.metagenes[o])] <- NA
      }
    }

    core[which(!is.na(core))] <- -1
    spot.i <- 0

    while (nrow(which(core == -1, arr.ind=T)) > 0)
    {
      start.pix <- which(core == -1, arr.ind=T)[1,]
      spot.i <- spot.i + 1
      core <- col.pix(core, start.pix[1], start.pix[2], spot.i, preferences$dim.1stLvlSom)
    }

    # define spot area around cores
    for (s.i in 1:max(core,na.rm=T))
    {
      n.sample.modules <- n.sample.modules + 1

      sample.spot.core.list[[n.sample.modules]] <-
        matrix(NA, preferences$dim.1stLvlSom, preferences$dim.1stLvlSom)

      sample.spot.core.list[[n.sample.modules]][which(core == s.i)] <-
        metadata[which(core == s.i),m]

      spot <- matrix(NA, preferences$dim.1stLvlSom, preferences$dim.1stLvlSom)
      spot[which(metadata.scaled[,m] > preferences$spot.threshold.modules)] <- -1

      start.pix <- which(!is.na(sample.spot.core.list[[n.sample.modules]]), arr.ind=T)
      start.pix <- start.pix[which.max(sample.spot.core.list[[n.sample.modules]][start.pix]),]

      spot <- col.pix(spot, start.pix[1], start.pix[2], 1, preferences$dim.1stLvlSom)

      sample.spot.list[[n.sample.modules]] <- matrix(NA, preferences$dim.1stLvlSom, preferences$dim.1stLvlSom)
      sample.spot.list[[n.sample.modules]][which(spot == 1)] <- metadata.scaled[which(spot == 1),m]
    }
  }

  # filter
  remove <- c()

  for (i in 1:n.sample.modules)
  {
    if (sum(!is.na(sample.spot.list[[i]])) <= 1)
    {
      # empty, i.e. does not exceed threshold -> remove
      remove <- c(remove, i)
    } else if (sum(!is.na(sample.spot.list[[i]])) < sum(!is.na(sample.spot.core.list[[i]])))
    {
      # core larger than spot -> shrink core
      sample.spot.core.list[[i]] <- sample.spot.list[[i]]
    }
  }

  if (length(remove) > 0)
  {
    sample.spot.list <- sample.spot.list[-remove]
    sample.spot.core.list <- sample.spot.core.list[-remove]
    n.sample.modules <- length(sample.spot.list)
  }

  o <- order(sapply(sample.spot.core.list, function(x) { mean(x, na.rm=T) }), decreasing=T)
  sample.spot.list <- sample.spot.list[o]
  sample.spot.core.list <- sample.spot.core.list[o]

  ## merge overlapping sample cores ##
  merged <- T

  while (merged)
  {
    merged <- F
    i <- 1

    while (i < length(sample.spot.list))
    {
      j <- i + 1

      while (j <= length(sample.spot.list))
      {
        core1 <- which(!is.na(sample.spot.core.list[[i]]))
        core2 <- which(!is.na(sample.spot.core.list[[j]]))

        if (any(core1 %in% core2))
        {
          # merge cores
          if (length(setdiff(core1,core2)) > 0)
          {
            sample.spot.core.list[[i]][setdiff(core2, core1)] <-
              sample.spot.core.list[[j]][setdiff(core2, core1)]
          }

          # merge spots
          spot1 <- which(!is.na(sample.spot.list[[i]]))
          spot2 <- which(!is.na(sample.spot.list[[j]]))

          if (length(setdiff(spot2,spot1)) > 0)
          {
            sample.spot.list[[i]][setdiff(spot2, spot1)] <-
              sample.spot.list[[j]][setdiff(spot2, spot1)]
          }

          # remove j
          sample.spot.list <- sample.spot.list[-j]
          sample.spot.core.list <- sample.spot.core.list[-j]

          merged <- T
        } else
        {
          j <- j + 1
        }
      }

      i <- i + 1
    }
  }

  o <- order(sapply(sample.spot.core.list, function(x) { mean(x,na.rm=T) }), decreasing=T)
  sample.spot.list <- sample.spot.list[o]
  sample.spot.core.list <- sample.spot.core.list[o]

  ## shrinking overlapping spots ##
  if (length(sample.spot.list) > 1)
  {
    for (i in 1:(length(sample.spot.list)-1))
    {
      for (j in (i+1):length(sample.spot.list))
      {
        spot1 <- which(!is.na(sample.spot.list[[i]]))
        spot2 <- which(!is.na(sample.spot.list[[j]]))

        if (any(spot1 %in% spot2))
        {
          spot12intersect <- which(spot1 %in% spot2)

          spot1 <- which(!is.na(sample.spot.list[[i]]), arr.ind=T)
          spot2 <- which(!is.na(sample.spot.list[[j]]), arr.ind=T)

          spot12intersect <- spot1[spot12intersect,,drop=F]

          core1.center <- colMeans(apply(which(!is.na(sample.spot.core.list[[i]]), arr.ind=T), 2, range))
          core2.center <- colMeans(apply(which(!is.na(sample.spot.core.list[[j]]), arr.ind=T), 2, range))

          spot12assoc <- apply(spot12intersect, 1, function(x)
          {
            which.min(c(sum((core1.center - x) ^ 2), sum((core2.center - x) ^ 2)))
          })

          sample.spot.list[[j]][spot12intersect[which(spot12assoc==1),,drop=F]] <- NA
          sample.spot.list[[i]][spot12intersect[which(spot12assoc==2),,drop=F]] <- NA
        }
      }
    }
  }

  ## define overexpression spots ##
  spot.list.overexpression <<- list()

  spot.list.overexpression$overview.map <<-
    apply(apply(metadata, 2, function(x) { (x - min(x)) / (max(x) - min(x)) }), 1, max)

  spot.list.overexpression$overview.mask <<- rep(NA, preferences$dim.1stLvlSom ^ 2)
  spot.list.overexpression$filtered <<- F
  spot.list.overexpression$spots <<- list()

  for (i in 1:length(sample.spot.list))
  {
    spot.list.overexpression$overview.mask[which(!is.na(sample.spot.list[[i]]))] <<- i

    spot.metagenes <- which(!is.na(sample.spot.list[[i]]))
    spot.genes <- rownames(indata)[which(som.nodes %in% spot.metagenes)]

    spot.list.overexpression$spots[[LETTERS[i]]] <<- list()
    spot.list.overexpression$spots[[LETTERS[i]]]$metagenes <<- spot.metagenes
    spot.list.overexpression$spots[[LETTERS[i]]]$genes <<- spot.genes
    spot.list.overexpression$spots[[LETTERS[i]]]$mask <<- rep(NA, preferences$dim.1stLvlSom * preferences$dim.1stLvlSom)
    spot.list.overexpression$spots[[LETTERS[i]]]$mask[spot.metagenes] <<- 1

    spot.list.overexpression$spots[[LETTERS[i]]]$position <<-
      colMeans(apply(som.result$code.sum[spot.metagenes, 1:2] + 1, 2, range))

    spot.list.overexpression$spots[[LETTERS[i]]]$beta.statistic <<-
      get.beta.statistic(set.data=metadata[spot.list.overexpression$spots[[LETTERS[i]]]$metagenes,,drop=F],
                         weights=som.result$code.sum[spot.list.overexpression$spots[[LETTERS[i]]]$metagenes,]$nobs)
  }

  o <- order(sapply(spot.list.overexpression$spots, function(x)
  {
    which.max(apply(metadata[x$metagenes,,drop=F], 2, mean))
  }))

  spot.list.overexpression$spots <<- spot.list.overexpression$spots[o]
  names(spot.list.overexpression$spots) <<- LETTERS[1:length(spot.list.overexpression$spots)]

  spot.list.overexpression$overview.mask[!is.na(spot.list.overexpression$overview.mask)] <<-
    match(spot.list.overexpression$overview.mask[!is.na(spot.list.overexpression$overview.mask)], o)

  spot.list.overexpression$spotdata <<-
    t(sapply(spot.list.overexpression$spots, function(x)
    {
      if (length(x$genes > 0))
      {
        colMeans(indata[x$genes,,drop=F])
      } else
      {
        rep(0, ncol(indata))
      }
    }))

  colnames(spot.list.overexpression$spotdata) <<- colnames(indata)


  ##### Underexpression Spots ######

  ## extract sample modules ##
  sample.spot.list <- list()
  sample.spot.core.list <- list()
  n.sample.modules <- 0

  for (m in 1:ncol(indata))
  {
    # define bigger core regions
    core <- matrix(NA, preferences$dim.1stLvlSom, preferences$dim.1stLvlSom)
    core[which(metadata.scaled[,m] < 1 - preferences$spot.threshold.modules)] <- -1

    spot.i <- 0

    while (nrow(which(core == -1, arr.ind=T)) > 0)
    {
      start.pix <- which(core == -1, arr.ind=T)[1,]
      spot.i <- spot.i + 1
      core <- col.pix(core, start.pix[1], start.pix[2], spot.i, preferences$dim.1stLvlSom)
    }

    # shrink each separate region to core size
    for (s.i in 1:max(core,na.rm=T))
    {
      if (sum(core == s.i, na.rm=T) > preferences$spot.coresize.modules)
      {
        core.metagenes <- which(core == s.i)
        o <- order(metadata[core.metagenes,m], decreasing=F)[1:preferences$spot.coresize.modules]
        core[setdiff(core.metagenes, core.metagenes[o])] <- NA
      }
    }

    core[which(!is.na(core))] <- -1
    spot.i <- 0

    while (nrow(which(core == -1, arr.ind=T)) > 0)
    {
      start.pix <- which(core == -1, arr.ind=T)[1,]
      spot.i <- spot.i + 1
      core <- col.pix(core, start.pix[1], start.pix[2], spot.i, preferences$dim.1stLvlSom)
    }

    # define spot area around cores
    for (s.i in 1:max(core,na.rm=T))
    {
      n.sample.modules <- n.sample.modules + 1
      sample.spot.core.list[[n.sample.modules]] <- matrix(NA, preferences$dim.1stLvlSom, preferences$dim.1stLvlSom)
      sample.spot.core.list[[n.sample.modules]][which(core == s.i)] <- metadata[which(core == s.i),m]

      spot <- matrix(NA, preferences$dim.1stLvlSom, preferences$dim.1stLvlSom)
      spot[which(metadata.scaled[,m] < 1 - preferences$spot.threshold.modules)] <- -1

      start.pix <- which(!is.na(sample.spot.core.list[[n.sample.modules]]), arr.ind=T)
      start.pix <- start.pix[which.max(sample.spot.core.list[[n.sample.modules]][start.pix]),]

      spot <- col.pix(spot, start.pix[1], start.pix[2], 1, preferences$dim.1stLvlSom)


      sample.spot.list[[n.sample.modules]] <- matrix(NA, preferences$dim.1stLvlSom, preferences$dim.1stLvlSom)
      sample.spot.list[[n.sample.modules]][which(spot == 1)] <- metadata.scaled[which(spot == 1),m]
    }
  }

  # filter
  remove <- c()

  for (i in 1:n.sample.modules)
  {
    if (sum(!is.na(sample.spot.list[[i]])) <= 1)
    {
      # empty, i.e. does not exceed threshold -> remove
      remove <- c(remove, i)
    } else if (sum(!is.na(sample.spot.list[[i]])) < sum(!is.na(sample.spot.core.list[[i]])))
    {
      # core larger than spot -> shrink core
      sample.spot.core.list[[i]] <- sample.spot.list[[i]]
    }
  }

  if (length(remove) > 0)
  {
    sample.spot.list <- sample.spot.list[-remove]
    sample.spot.core.list <- sample.spot.core.list[-remove]
    n.sample.modules <- length(sample.spot.list)
  }

  o <- order(sapply(sample.spot.core.list, function(x) { mean(x,na.rm=T) }), decreasing=F)
  sample.spot.list <- sample.spot.list[o]
  sample.spot.core.list <- sample.spot.core.list[o]

  ## merge overlapping sample cores ##
  merged <- T

  while (merged)
  {
    merged <- F
    i <- 1

    while (i < length(sample.spot.list))
    {
      j <- i + 1

      while (j <= length(sample.spot.list))
      {
        core1 <- which(!is.na(sample.spot.core.list[[i]]))
        core2 <- which(!is.na(sample.spot.core.list[[j]]))

        if (any(core1 %in% core2))
        {
          # merge cores
          if (length(setdiff(core1,core2)) > 0)
          {
            sample.spot.core.list[[i]][setdiff(core2, core1)] <-
              sample.spot.core.list[[j]][setdiff(core2,core1)]
          }

          # merge spots
          spot1 <- which(!is.na(sample.spot.list[[i]]))
          spot2 <- which(!is.na(sample.spot.list[[j]]))

          if (length(setdiff(spot2,spot1)) > 0)
          {
            sample.spot.list[[i]][setdiff(spot2, spot1)] <-
              sample.spot.list[[j]][setdiff(spot2, spot1)]
          }

          # remove j
          sample.spot.list <- sample.spot.list[-j]
          sample.spot.core.list <- sample.spot.core.list[-j]

          merged <- T
        } else
        {
          j <- j + 1
        }
      }

      i <- i + 1
    }

  }

  o <- order(sapply(sample.spot.core.list,function(x) { mean(x,na.rm=T) }), decreasing=F)
  sample.spot.list <- sample.spot.list[o]
  sample.spot.core.list <- sample.spot.core.list[o]

  ## shrinking overlapping spots ##
  if (length(sample.spot.list) > 1)
  {
    for (i in 1:(length(sample.spot.list)-1))
    {
      for (j in (i+1):length(sample.spot.list))
      {
        spot1 <- which(!is.na(sample.spot.list[[i]]))
        spot2 <- which(!is.na(sample.spot.list[[j]]))


        if (any(spot1 %in% spot2))
        {
          spot12intersect <- which(spot1 %in% spot2)

          spot1 <- which(!is.na(sample.spot.list[[i]]), arr.ind=T)
          spot2 <- which(!is.na(sample.spot.list[[j]]), arr.ind=T)

          spot12intersect <- spot1[spot12intersect,,drop=F]

          core1.center <- colMeans(apply(which(!is.na(sample.spot.core.list[[i]]), arr.ind=T), 2, range))
          core2.center <- colMeans(apply(which(!is.na(sample.spot.core.list[[j]]), arr.ind=T), 2, range))

          spot12assoc <- apply(spot12intersect, 1, function(x)
          {
            which.min(c(sum((core1.center - x) ^ 2), sum((core2.center - x) ^ 2)))
          })

          sample.spot.list[[j]][spot12intersect[which(spot12assoc==1),,drop=F]] <- NA
          sample.spot.list[[i]][spot12intersect[which(spot12assoc==2),,drop=F]] <- NA
        }
      }
    }
  }

  ## define underexpression spots ##
  spot.list.underexpression <<- list()

  spot.list.underexpression$overview.map <<-
    apply(apply(metadata, 2, function(x) { (x - min(x)) / (max(x) - min(x)) }), 1, min)

  spot.list.underexpression$overview.mask <<- rep(NA, preferences$dim.1stLvlSom ^ 2)
  spot.list.underexpression$filtered <<- F
  spot.list.underexpression$spots <<- list()

  for (i in 1:length(sample.spot.list))
  {
    spot.list.underexpression$overview.mask[which(!is.na(sample.spot.list[[i]]))] <<- i

    spot.metagenes <- which(!is.na(sample.spot.list[[i]]))
    spot.genes <- rownames(indata)[which(som.nodes %in% spot.metagenes)]

    spot.list.underexpression$spots[[letters[i]]] <<- list()
    spot.list.underexpression$spots[[letters[i]]]$metagenes <<- spot.metagenes
    spot.list.underexpression$spots[[letters[i]]]$genes <<- spot.genes
    spot.list.underexpression$spots[[letters[i]]]$mask <<- rep(NA, preferences$dim.1stLvlSom * preferences$dim.1stLvlSom)
    spot.list.underexpression$spots[[letters[i]]]$mask[spot.metagenes] <<- 1

    spot.list.underexpression$spots[[letters[i]]]$position <<-
      colMeans(apply(som.result$code.sum[spot.metagenes, 1:2]+1, 2, range))

    spot.list.underexpression$spots[[letters[i]]]$beta.statistic <<-
      get.beta.statistic(set.data=metadata[spot.list.underexpression$spots[[letters[i]]]$metagenes,,drop=F],
                         weights=som.result$code.sum[spot.list.underexpression$spots[[letters[i]]]$metagenes,]$nobs)

  }

  o <- order(sapply(spot.list.underexpression$spots, function(x)
  {
    which.min(apply(metadata[x$metagenes,,drop=F], 2, mean))
  }))

  spot.list.underexpression$spots <<- spot.list.underexpression$spots[o]
  names(spot.list.underexpression$spots) <<- letters[1:length(spot.list.underexpression$spots)]

  spot.list.underexpression$overview.mask[!is.na(spot.list.underexpression$overview.mask)] <<-
    match(spot.list.underexpression$overview.mask[!is.na(spot.list.underexpression$overview.mask)], o)

  spot.list.underexpression$spotdata <<-
    t(sapply(spot.list.underexpression$spots, function(x)
    {
      if (length(x$genes > 0))
      {
        colMeans(indata[x$genes,,drop=F])
      } else
      {
        rep(0, ncol(indata))
      }
    }))

  colnames(spot.list.underexpression$spotdata) <<- colnames(indata)


  ##### Correlation Cluster ######
  spot.list.correlation <<- list()
  spot.list.correlation$overview.map <<- NA
  spot.list.correlation$overview.mask <<- rep(NA, preferences$dim.1stLvlSom ^ 2)
  spot.list.correlation$filtered <<- F
  spot.list.correlation$spots <<- list()

  c.map <- cor(t(metadata))
  diag(c.map) <- NA
  rownames(c.map) <- c(1:(preferences$dim.1stLvlSom * preferences$dim.1stLvlSom))
  colnames(c.map) <- c(1:(preferences$dim.1stLvlSom * preferences$dim.1stLvlSom))

  c.cluster <- rep(NA, preferences$dim.1stLvlSom * preferences$dim.1stLvlSom)
  names(c.cluster) <- c(1:(preferences$dim.1stLvlSom * preferences$dim.1stLvlSom))

  count.cluster <- 1

  while (is.matrix(c.map) && nrow(c.map) > 0 && ncol(c.map) > 0 && count.cluster <= 30)
  {
    start.node <- rownames(c.map)[which(c.map == max(c.map, na.rm=T), arr.ind=T)[1,1]]

    cluster <- names(which(c.map[start.node,] > 0.90))
    cluster <- c(start.node, cluster)

    if (length(cluster) >= preferences$dim.1stLvlSom / 2)
    {
      c.cluster[cluster] <- count.cluster
      geneset.genes <- rownames(indata)[which(som.nodes %in% as.numeric(cluster))]

      spot.list.correlation$overview.mask[as.numeric(cluster)] <<- count.cluster
      spot.list.correlation$spots[[count.cluster]] <<- list()
      spot.list.correlation$spots[[count.cluster]]$metagenes <<- as.numeric(cluster)
      spot.list.correlation$spots[[count.cluster]]$genes <<- geneset.genes
      spot.list.correlation$spots[[count.cluster]]$mask <<- rep(NA, preferences$dim.1stLvlSom * preferences$dim.1stLvlSom)
      spot.list.correlation$spots[[count.cluster]]$mask[as.numeric(cluster)] <<- 1

      spot.list.correlation$spots[[count.cluster]]$position <<-
        apply(apply(som.result$code.sum[cluster, 1:2], 2, range), 2, mean) + 0.5

      spot.list.correlation$spots[[count.cluster]]$beta.statistic <<-
        get.beta.statistic(set.data=metadata[spot.list.correlation$spots[[count.cluster]]$metagenes,,drop=F],
                           weights=som.result$code.sum[spot.list.correlation$spots[[count.cluster]]$metagenes,]$nobs)

      count.cluster <- count.cluster + 1
    }

    c.map <- c.map[-which(rownames(c.map) %in% cluster), -which(colnames(c.map) %in% cluster)]
  }

  spot.list.correlation$overview.map <<- c.cluster

  o <- order(sapply(spot.list.correlation$spots, function(x)
  {
    which.max(apply(metadata[x$metagenes,,drop=F], 2, mean))
  }))

  spot.list.correlation$spots <<- spot.list.correlation$spots[o]
  names(spot.list.correlation$spots) <<- LETTERS[1:length(spot.list.correlation$spots)]

  spot.list.correlation$overview.mask[!is.na(spot.list.correlation$overview.mask)] <<-
    match(spot.list.correlation$overview.mask[!is.na(spot.list.correlation$overview.mask)], o)

  spot.list.correlation$spotdata <<-
    t(sapply(spot.list.correlation$spots, function(x)
    {
      if (length(x$genes > 0))
      {
        colMeans(indata[x$genes,,drop=F])
      } else
      {
        rep(0, ncol(indata))
      }
    }))

  colnames(spot.list.correlation$spotdata) <<- colnames(indata)


  ##### K-Means Clustering #####
  n.cluster <- preferences$dim.1stLvlSom / 2
  prototypes <- metadata[round(seq(1, preferences$dim.1stLvlSom^2, length.out=n.cluster)),]
  res <- kmeans(metadata, prototypes)

  spot.list.kmeans <<- list()
  spot.list.kmeans$overview.map <<- matrix(res$cluster, preferences$dim.1stLvlSom, preferences$dim.1stLvlSom)
  spot.list.kmeans$overview.mask <<- rep(NA, preferences$dim.1stLvlSom ^ 2)
  spot.list.kmeans$filtered <<- F
  spot.list.kmeans$spots <<- list()

  for (i in 1:n.cluster)
  {
    nodes <- which(res$cluster == i)
    geneset.genes <- rownames(indata)[which(som.nodes %in% nodes)]

    spot.list.kmeans$overview.mask[as.numeric(nodes)] <<- i
    spot.list.kmeans$spots[[LETTERS[i]]] <<- list()
    spot.list.kmeans$spots[[LETTERS[i]]]$metagenes <<- as.numeric(nodes)
    spot.list.kmeans$spots[[LETTERS[i]]]$genes <<- geneset.genes
    spot.list.kmeans$spots[[LETTERS[i]]]$mask <<- rep(NA, preferences$dim.1stLvlSom * preferences$dim.1stLvlSom)
    spot.list.kmeans$spots[[LETTERS[i]]]$mask[as.numeric(nodes)] <<- 1

    spot.list.kmeans$spots[[LETTERS[i]]]$position <<-
      apply(apply(som.result$code.sum[nodes, 1:2], 2, range), 2, mean) + 0.5

    spot.list.kmeans$spots[[LETTERS[i]]]$beta.statistic <<-
      get.beta.statistic(set.data=metadata[spot.list.kmeans$spots[[LETTERS[i]]]$metagenes,,drop=F],
                         weights=som.result$code.sum[spot.list.kmeans$spots[[LETTERS[i]]]$metagenes,]$nobs)
  }

  o <- order(sapply(spot.list.kmeans$spots, function(x)
  {
    which.max(apply(metadata[x$metagenes,,drop=F], 2, mean))
  }))

  spot.list.kmeans$spots <<- spot.list.kmeans$spots[o]
  names(spot.list.kmeans$spots) <<- LETTERS[1:length(spot.list.kmeans$spots)]

  spot.list.kmeans$overview.mask[!is.na(spot.list.kmeans$overview.mask)] <<-
    match(spot.list.kmeans$overview.mask[!is.na(spot.list.kmeans$overview.mask)], o)

  spot.list.kmeans$spotdata <<-
    t(sapply(spot.list.kmeans$spots, function(x)
    {
      if (length(x$genes > 0))
      {
        colMeans(indata[x$genes,,drop=F])
      } else
      {
        rep(0, ncol(indata))
      }
    }))

  colnames(spot.list.kmeans$spotdata) <<- colnames(indata)


  ##### Group Spots ######
  if (length(unique(group.labels)) > 1)
  {
    group.metadata.scaled <- apply(group.metadata, 2, function(x)
    {
      (x-min(x)) / (max(x)-min(x))
    })

    ## extract sample modules ##
    sample.spot.list <- list()
    sample.spot.core.list <- list()
    n.sample.modules <- 0

    for (m in 1:ncol(group.metadata))
    {
      # define bigger core regions
      core <- matrix(NA, preferences$dim.1stLvlSom, preferences$dim.1stLvlSom)
      core[which(group.metadata.scaled[,m] > preferences$spot.threshold.groupmap)] <- -1

      spot.i <- 0

      while (nrow(which(core == -1, arr.ind=T)) > 0)
      {
        start.pix <- which(core == -1, arr.ind=T)[1,]

        spot.i <- spot.i + 1
        core <- col.pix(core, start.pix[1], start.pix[2], spot.i, preferences$dim.1stLvlSom)
      }

      # shrink each separate region to core size
      for (s.i in 1:max(core,na.rm=T))
      {
        if (sum(core == s.i, na.rm=T) > preferences$spot.coresize.groupmap)
        {
          core.metagenes <- which(core == s.i)
          o <- order(group.metadata[core.metagenes,m], decreasing=T)[1:preferences$spot.coresize.groupmap]
          core[setdiff(core.metagenes, core.metagenes[o])] <- NA
        }
      }

      core[which(!is.na(core))] <- -1
      spot.i <- 0

      while (nrow(which(core == -1, arr.ind=T)) > 0)
      {
        start.pix <- which(core == -1, arr.ind=T)[1,]
        spot.i <- spot.i + 1
        core <- col.pix(core, start.pix[1], start.pix[2], spot.i, preferences$dim.1stLvlSom)
      }

      # define spot area around cores
      for (s.i in 1:max(core,na.rm=T))
      {
        n.sample.modules <- n.sample.modules + 1

        sample.spot.core.list[[n.sample.modules]] <- matrix(NA, preferences$dim.1stLvlSom, preferences$dim.1stLvlSom)
        sample.spot.core.list[[n.sample.modules]][which(core == s.i)] <- metadata[which(core == s.i),m]

        spot <- matrix(NA, preferences$dim.1stLvlSom, preferences$dim.1stLvlSom)
        spot[which(group.metadata.scaled[,m] > preferences$spot.threshold.groupmap)] <- -1

        start.pix <- which(!is.na(sample.spot.core.list[[n.sample.modules]]), arr.ind=T)
        start.pix <- start.pix[which.max(sample.spot.core.list[[n.sample.modules]][start.pix]),]

        spot <- col.pix(spot, start.pix[1], start.pix[2], 1, preferences$dim.1stLvlSom)

        sample.spot.list[[n.sample.modules]] <- matrix(NA, preferences$dim.1stLvlSom, preferences$dim.1stLvlSom)
        sample.spot.list[[n.sample.modules]][which(spot == 1)] <- group.metadata.scaled[which(spot == 1),m]
      }
    }

    # filter
    remove <- c()

    for (i in 1:n.sample.modules)
    {
      if (sum(!is.na(sample.spot.list[[i]])) <= 1)
      {
        # empty, i.e. does not exceed threshold -> remove
        remove <- c(remove, i)
      } else if (sum(!is.na(sample.spot.list[[i]])) < sum(!is.na(sample.spot.core.list[[i]])))
      {
        # core larger than spot -> shrink core
        sample.spot.core.list[[i]] <- sample.spot.list[[i]]
      }
    }

    if (length(remove) > 0)
    {
      sample.spot.list <- sample.spot.list[-remove]
      sample.spot.core.list <- sample.spot.core.list[-remove]
      n.sample.modules <- length(sample.spot.list)
    }

    o <- order(sapply(sample.spot.core.list,function(x) mean(x,na.rm=T)), decreasing=T)
    sample.spot.list <- sample.spot.list[o]
    sample.spot.core.list <- sample.spot.core.list[o]

    ## merge overlapping sample cores ##
    merged <- T

    while (merged)
    {
      merged <- F
      i <- 1

      while (i < length(sample.spot.list))
      {
        j <- i + 1

        while (j <= length(sample.spot.list))
        {
          core1 <- which(!is.na(sample.spot.core.list[[i]]))
          core2 <- which(!is.na(sample.spot.core.list[[j]]))

          if (any(core1 %in% core2))
          {
            # merge cores
            if (length(setdiff(core1,core2)) > 0)
            {
              sample.spot.core.list[[i]][setdiff(core2,core1)] <-
                sample.spot.core.list[[j]][setdiff(core2,core1)]
            }

            # merge spots
            spot1 <- which(!is.na(sample.spot.list[[i]]))
            spot2 <- which(!is.na(sample.spot.list[[j]]))

            if (length(setdiff(spot2,spot1)) > 0)
            {
              sample.spot.list[[i]][setdiff(spot2,spot1)] <-
                sample.spot.list[[j]][setdiff(spot2,spot1)]
            }

            # remove j
            sample.spot.list <- sample.spot.list[-j]
            sample.spot.core.list <- sample.spot.core.list[-j]

            merged <- T
          } else
          {
            j <- j + 1
          }
        }

        i <- i + 1
      }
    }

    o <- order(sapply(sample.spot.core.list, function(x) { mean(x, na.rm=T) }), decreasing=T)
    sample.spot.list <- sample.spot.list[o]
    sample.spot.core.list <- sample.spot.core.list[o]

    ## shrinking overlapping spots ##
    if (length(sample.spot.list) > 1)
    {
      for (i in 1:(length(sample.spot.list)-1))
      {
        for (j in (i+1):length(sample.spot.list))
        {
          spot1 <- which(!is.na(sample.spot.list[[i]]))
          spot2 <- which(!is.na(sample.spot.list[[j]]))

          if (any(spot1 %in% spot2))
          {
            spot12intersect <- which(spot1 %in% spot2)

            spot1 <- which(!is.na(sample.spot.list[[i]]), arr.ind=T)
            spot2 <- which(!is.na(sample.spot.list[[j]]), arr.ind=T)

            spot12intersect <- spot1[spot12intersect, ,drop=F]

            core1.center <- colMeans(apply(which(!is.na(sample.spot.core.list[[i]]), arr.ind=T), 2, range))
            core2.center <- colMeans(apply(which(!is.na(sample.spot.core.list[[j]]), arr.ind=T), 2, range))

            spot12assoc <- apply(spot12intersect, 1, function(x)
            {
              which.min(c(sum((core1.center - x) ^ 2), sum((core2.center - x) ^ 2)))
            })

            sample.spot.list[[j]][spot12intersect[which(spot12assoc==1),,drop=F]] <- NA
            sample.spot.list[[i]][spot12intersect[which(spot12assoc==2),,drop=F]] <- NA
          }
        }
      }
    }

    ## define overexpression spots ##
    spot.list.group.overexpression <<- list()

    spot.list.group.overexpression$overview.map <<-
      apply(apply(group.metadata, 2, function(x){ (x - min(x)) / (max(x) - min(x)) }), 1, max)

    spot.list.group.overexpression$overview.mask <<- rep(NA, preferences$dim.1stLvlSom ^ 2)
    spot.list.group.overexpression$filtered <<- F
    spot.list.group.overexpression$spots <<- list()

    for (i in 1:length(sample.spot.list))
    {
      spot.list.group.overexpression$overview.mask[which(!is.na(sample.spot.list[[i]]))] <<- i

      spot.metagenes <- which(!is.na(sample.spot.list[[i]]))
      spot.genes <- rownames(indata)[which(som.nodes %in% spot.metagenes)]

      spot.list.group.overexpression$spots[[LETTERS[i]]] <<- list()
      spot.list.group.overexpression$spots[[LETTERS[i]]]$metagenes <<- spot.metagenes
      spot.list.group.overexpression$spots[[LETTERS[i]]]$genes <<- spot.genes
      spot.list.group.overexpression$spots[[LETTERS[i]]]$mask <<- rep(NA, preferences$dim.1stLvlSom * preferences$dim.1stLvlSom)
      spot.list.group.overexpression$spots[[LETTERS[i]]]$mask[spot.metagenes] <<- 1

      spot.list.group.overexpression$spots[[LETTERS[i]]]$position <<-
        colMeans(apply(som.result$code.sum[spot.metagenes, 1:2]+1, 2, range))

      spot.list.group.overexpression$spots[[LETTERS[i]]]$beta.statistic <<-
        get.beta.statistic(set.data=metadata[spot.list.group.overexpression$spots[[LETTERS[i]]]$metagenes,,drop=F],
                           weights=som.result$code.sum[spot.list.group.overexpression$spots[[LETTERS[i]]]$metagenes,]$nobs)
    }

    start.spot <- which.min(sapply(spot.list.group.overexpression$spots, function(x)
    {
      which.max(apply(group.metadata[x$metagenes,,drop=F], 2, mean))
    }))

    spot.arcs <- sapply(spot.list.group.overexpression$spots, function(x)
    {
      -atan2(x$position['y'] - preferences$dim.1stLvlSom / 2, x$position['x'] - preferences$dim.1stLvlSom / 2)
    })

    spot.arcs <- spot.arcs - spot.arcs[start.spot]

    if (any(spot.arcs<0))
    {
      spot.arcs[which(spot.arcs<0)] <- spot.arcs[which(spot.arcs<0)] + (2 * pi)
    }

    o <- order(spot.arcs)

    spot.list.group.overexpression$spots <<- spot.list.group.overexpression$spots[o]
    names(spot.list.group.overexpression$spots) <<- LETTERS[1:length(spot.list.group.overexpression$spots)]

    spot.list.group.overexpression$overview.mask[!is.na(spot.list.group.overexpression$overview.mask)] <<-
      match(spot.list.group.overexpression$overview.mask[!is.na(spot.list.group.overexpression$overview.mask)], o)

    spot.list.group.overexpression$spotdata <<-
      t(sapply(spot.list.group.overexpression$spots, function(x)
      {
        if (length(x$genes > 0))
        {
          colMeans(indata[x$genes,,drop=F])
        } else
        {
          rep(0, ncol(indata))
        }
      }))

    colnames(spot.list.group.overexpression$spotdata) <<- colnames(indata)
  }
}
