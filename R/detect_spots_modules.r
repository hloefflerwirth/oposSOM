pipeline.detectSpotsModules <- function(env)
{
  env <- pipeline.detectOverexpressionModules(env)
  env <- pipeline.detectUnderexpressionModules(env)
  env <- pipeline.detectCorrelationModules(env)
  env <- pipeline.detectKMeansModules(env)
  env <- pipeline.detectGroupOverexpressionModules(env)
  env <- pipeline.detectDMapModules(env)
  
  # check standard spot modules
  if( length( env[[paste("spot.list.", env$preferences$standard.spot.modules,sep="")]]$spots ) < 2 )
  {
    env$preferences$standard.spot.modules <- "kmeans"
    util.warn("Invalid value of \"standard.spot.modules\": Too few spots detected. Using \"kmeans\"")
  } 
  
  return(env)
}
  
  
pipeline.detectOverexpressionModules <- function(env)
{
  metadata.scaled <- apply(env$metadata, 2, function(x)
  {
    (x-min(x)) / (max(x)-min(x))
  })
  
  # extract sample modules
  sample.spot.list <- list()
  sample.spot.core.list <- list()
  n.sample.modules <- 0

  for (m in 1:ncol(env$indata))
  {
    # define bigger core regions
    core <- matrix(NA, env$preferences$dim.1stLvlSom, env$preferences$dim.1stLvlSom)
    core[which(metadata.scaled[,m] > env$preferences$spot.threshold.modules)] <- -1

    spot.i = 0

    while (nrow(which(core == -1, arr.ind=TRUE)) > 0)
    {
      start.pix <- which(core == -1, arr.ind=TRUE)[1,]
      spot.i <- spot.i + 1
      core <- col.pix(core, start.pix[1], start.pix[2], spot.i, env$preferences$dim.1stLvlSom)
    }
    
    # shrink each separate region to core size
    for (s.i in 1:max(core, na.rm=TRUE))
    {
      if (sum(core == s.i, na.rm=TRUE) > env$preferences$spot.coresize.modules)
      {
        core.metagenes <- which(core == s.i)
        o <- order(env$metadata[core.metagenes,m], decreasing=TRUE)[1:env$preferences$spot.coresize.modules]
        core[setdiff(core.metagenes, core.metagenes[o])] <- NA
      }
    }

    core[which(!is.na(core))] <- -1
    spot.i <- 0

    while (nrow(which(core == -1, arr.ind=TRUE)) > 0)
    {
      start.pix <- which(core == -1, arr.ind=TRUE)[1,]
      spot.i <- spot.i + 1
      core <- col.pix(core, start.pix[1], start.pix[2], spot.i, env$preferences$dim.1stLvlSom)
    }
    
    # define spot area around cores
    for (s.i in 1:max(core,na.rm=TRUE))
    {
      n.sample.modules <- n.sample.modules + 1

      sample.spot.core.list[[n.sample.modules]] <-
        matrix(NA, env$preferences$dim.1stLvlSom, env$preferences$dim.1stLvlSom)

      sample.spot.core.list[[n.sample.modules]][which(core == s.i)] <-
        env$metadata[which(core == s.i),m]

      spot <- matrix(NA, env$preferences$dim.1stLvlSom, env$preferences$dim.1stLvlSom)
      spot[which(metadata.scaled[,m] > env$preferences$spot.threshold.modules)] <- -1

      start.pix <- which(!is.na(sample.spot.core.list[[n.sample.modules]]), arr.ind=TRUE)
      start.pix <- start.pix[which.max(sample.spot.core.list[[n.sample.modules]][start.pix]),]

      spot <- col.pix(spot, start.pix[1], start.pix[2], 1, env$preferences$dim.1stLvlSom)

      sample.spot.list[[n.sample.modules]] <- matrix(NA, env$preferences$dim.1stLvlSom, env$preferences$dim.1stLvlSom)
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

  o <- order(sapply(sample.spot.core.list, function(x) { mean(x, na.rm=TRUE) }), decreasing=TRUE)
  sample.spot.list <- sample.spot.list[o]
  sample.spot.core.list <- sample.spot.core.list[o]

  ## merge overlapping sample cores ##
  merged <- TRUE

  while (merged)
  {
    merged <- FALSE
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

          merged <- TRUE
        } else
        {
          j <- j + 1
        }
      }

      i <- i + 1
    }
  }
  
  o <- order(sapply(sample.spot.core.list, function(x) { mean(x,na.rm=TRUE) }), decreasing=TRUE)
  sample.spot.list <- sample.spot.list[o]
  sample.spot.core.list <- sample.spot.core.list[o]

  ## shrinking overlapping spots ##
  if (length(sample.spot.list) > 1)
  {
    for (i in seq_len(length(sample.spot.list)-1))
    {
      for (j in seq(i+1, length(sample.spot.list)))
      {
        spot1 <- which(!is.na(sample.spot.list[[i]]))
        spot2 <- which(!is.na(sample.spot.list[[j]]))

        if (any(spot1 %in% spot2))
        {
          spot12intersect <- which(spot1 %in% spot2)

          spot1 <- which(!is.na(sample.spot.list[[i]]), arr.ind=TRUE)
          spot2 <- which(!is.na(sample.spot.list[[j]]), arr.ind=TRUE)

          spot12intersect <- spot1[spot12intersect,,drop=FALSE]

          core1.center <- colMeans(apply(which(!is.na(sample.spot.core.list[[i]]), arr.ind=TRUE), 2, range))
          core2.center <- colMeans(apply(which(!is.na(sample.spot.core.list[[j]]), arr.ind=TRUE), 2, range))

          spot12assoc <- apply(spot12intersect, 1, function(x)
          {
            which.min(c(sum((core1.center - x) ^ 2), sum((core2.center - x) ^ 2)))
          })

          sample.spot.list[[j]][spot12intersect[which(spot12assoc==1),,drop=FALSE]] <- NA
          sample.spot.list[[i]][spot12intersect[which(spot12assoc==2),,drop=FALSE]] <- NA
        }
      }
    }
  }
  
  ## no spot detected? ##
  if( length(sample.spot.list)==0 )
  {
	  util.warn("No overexpression spot detectable.")
    sample.spot.list = list( matrix(NA, env$preferences$dim.1stLvlSom, env$preferences$dim.1stLvlSom) )
    sample.spot.list[[1]][ which.max(rowMeans(env$metadata)) ] = 1
  }

  ## check for split modules

  for (i in seq_along(sample.spot.list))
  {
    core <- sample.spot.list[[i]]
    core[which(!is.na(core))] <- -1
    spot.i <- 0
    
    while (nrow(which(core == -1, arr.ind=TRUE)) > 0)
    {
      start.pix <- which(core == -1, arr.ind=TRUE)[1,]
      spot.i <- spot.i + 1
      core <- col.pix(core, start.pix[1], start.pix[2], spot.i, env$preferences$dim.1stLvlSom)
    }
    
    if( max(core,na.rm=TRUE) > 1 )
    {
      sample.spot.list <- sample.spot.list[-i]
      
      for( module.i in 1:max(core,na.rm=TRUE) )
      {
        if( sum(core == module.i, na.rm=TRUE) >= env$preferences$spot.coresize.modules ) 
        {
          split.module <- core == module.i
          split.module[which(!split.module)] <- NA
          sample.spot.list[[length(sample.spot.list)+1]] <- split.module
        }
      }
    }
  }
  
  ## define overexpression spots ##
  env$spot.list.overexpression <- list()

  env$spot.list.overexpression$overview.map <-
    apply(apply(env$metadata, 2, function(x) { (x - min(x)) / (max(x) - min(x)) }), 1, max)

  env$spot.list.overexpression$overview.mask <- rep(NA, env$preferences$dim.1stLvlSom ^ 2)
  env$spot.list.overexpression$spots <- list()

  for (i in seq_along(sample.spot.list))
  {
    spot.metagenes <- which(!is.na(sample.spot.list[[i]]))
    spot.genes <- rownames(env$indata)[which(env$som.result$feature.BMU %in% spot.metagenes)]

    if (length(spot.genes) > 0)
    {
      env$spot.list.overexpression$overview.mask[which(!is.na(sample.spot.list[[i]]))] <- i
      env$spot.list.overexpression$spots[[env$LETTERS[i]]] <- list()
      env$spot.list.overexpression$spots[[env$LETTERS[i]]]$metagenes <- spot.metagenes
      env$spot.list.overexpression$spots[[env$LETTERS[i]]]$genes <- spot.genes
      env$spot.list.overexpression$spots[[env$LETTERS[i]]]$mask <- rep(NA, env$preferences$dim.1stLvlSom * env$preferences$dim.1stLvlSom)
      env$spot.list.overexpression$spots[[env$LETTERS[i]]]$mask[spot.metagenes] <- 1

      env$spot.list.overexpression$spots[[env$LETTERS[i]]]$position <-
        colMeans(apply(env$som.result$node.summary[spot.metagenes, 1:2] + 1, 2, range))

      env$spot.list.overexpression$spots[[env$LETTERS[i]]]$beta.statistic <-
        get.beta.statistic(set.data=env$metadata[env$spot.list.overexpression$spots[[env$LETTERS[i]]]$metagenes,,drop=FALSE],
                           weights=env$som.result$node.summary[env$spot.list.overexpression$spots[[env$LETTERS[i]]]$metagenes,]$n.features)
    }
  }

  o <- order(sapply(env$spot.list.overexpression$spots, function(x)
  {
    which.max(apply(env$metadata[x$metagenes,,drop=FALSE], 2, mean))
  }))
  
  env$spot.list.overexpression$spots <- env$spot.list.overexpression$spots[o]
  names(env$spot.list.overexpression$spots) <- env$LETTERS[seq_along(env$spot.list.overexpression$spots)]

  env$spot.list.overexpression$overview.mask[!is.na(env$spot.list.overexpression$overview.mask)] <-
    match(env$spot.list.overexpression$overview.mask[!is.na(env$spot.list.overexpression$overview.mask)], o)

  env$spot.list.overexpression$spotdata <-
    t(sapply(env$spot.list.overexpression$spots, function(x)
    {
      if (length(x$genes > 0))
      {
        colMeans(env$indata[x$genes,,drop=FALSE],na.rm=TRUE)
      } else
      {
        rep(0, ncol(env$indata))
      }
    }))

  colnames(env$spot.list.overexpression$spotdata) <- colnames(env$indata)

  return(env)   
}


pipeline.detectUnderexpressionModules <- function(env)
{
  metadata.scaled <- apply(env$metadata, 2, function(x)
  {
    (x-min(x)) / (max(x)-min(x))
  })
  
  ## extract sample modules ##
  sample.spot.list <- list()
  sample.spot.core.list <- list()
  n.sample.modules <- 0

  for (m in 1:ncol(env$indata))
  {
    # define bigger core regions
    core <- matrix(NA, env$preferences$dim.1stLvlSom, env$preferences$dim.1stLvlSom)
    core[which(metadata.scaled[,m] < 1 - env$preferences$spot.threshold.modules)] <- -1
    
    spot.i <- 0
    while (nrow(which(core == -1, arr.ind=TRUE)) > 0)
    {
      start.pix <- which(core == -1, arr.ind=TRUE)[1,]
      spot.i <- spot.i + 1
      core <- col.pix(core, start.pix[1], start.pix[2], spot.i, env$preferences$dim.1stLvlSom)
    }
    
    # shrink each separate region to core size
    for (s.i in 1:max(core,na.rm=TRUE))
    {
      if (sum(core == s.i, na.rm=TRUE) > env$preferences$spot.coresize.modules)
      {
        core.metagenes <- which(core == s.i)
        o <- order(env$metadata[core.metagenes,m], decreasing=FALSE)[1:env$preferences$spot.coresize.modules]
        core[setdiff(core.metagenes, core.metagenes[o])] <- NA
      }
    }

    core[which(!is.na(core))] <- -1
    spot.i <- 0

    while (nrow(which(core == -1, arr.ind=TRUE)) > 0)
    {
      start.pix <- which(core == -1, arr.ind=TRUE)[1,]
      spot.i <- spot.i + 1
      core <- col.pix(core, start.pix[1], start.pix[2], spot.i, env$preferences$dim.1stLvlSom)
    }
    
    # define spot area around cores
    for (s.i in 1:max(core,na.rm=TRUE))
    {
      n.sample.modules <- n.sample.modules + 1
      sample.spot.core.list[[n.sample.modules]] <- matrix(NA, env$preferences$dim.1stLvlSom, env$preferences$dim.1stLvlSom)
      sample.spot.core.list[[n.sample.modules]][which(core == s.i)] <- env$metadata[which(core == s.i),m]

      spot <- matrix(NA, env$preferences$dim.1stLvlSom, env$preferences$dim.1stLvlSom)
      spot[which(metadata.scaled[,m] < 1 - env$preferences$spot.threshold.modules)] <- -1

      start.pix <- which(!is.na(sample.spot.core.list[[n.sample.modules]]), arr.ind=TRUE)
      start.pix <- start.pix[which.max(sample.spot.core.list[[n.sample.modules]][start.pix]),]

      spot <- col.pix(spot, start.pix[1], start.pix[2], 1, env$preferences$dim.1stLvlSom)


      sample.spot.list[[n.sample.modules]] <- matrix(NA, env$preferences$dim.1stLvlSom, env$preferences$dim.1stLvlSom)
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

  o <- order(sapply(sample.spot.core.list, function(x) { mean(x,na.rm=TRUE) }), decreasing=FALSE)
  sample.spot.list <- sample.spot.list[o]
  sample.spot.core.list <- sample.spot.core.list[o]

  ## merge overlapping sample cores ##
  merged <- TRUE

  while (merged)
  {
    merged <- FALSE
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

          merged <- TRUE
        } else
        {
          j <- j + 1
        }
      }

      i <- i + 1
    }

  }

  o <- order(sapply(sample.spot.core.list,function(x) { mean(x,na.rm=TRUE) }), decreasing=FALSE)
  sample.spot.list <- sample.spot.list[o]
  sample.spot.core.list <- sample.spot.core.list[o]
  
  ## shrinking overlapping spots ##
  if (length(sample.spot.list) > 1)
  {
    for (i in seq_len(length(sample.spot.list)-1))
    {
      for (j in seq(i+1, length(sample.spot.list)))
      {
        spot1 <- which(!is.na(sample.spot.list[[i]]))
        spot2 <- which(!is.na(sample.spot.list[[j]]))


        if (any(spot1 %in% spot2))
        {
          spot12intersect <- which(spot1 %in% spot2)

          spot1 <- which(!is.na(sample.spot.list[[i]]), arr.ind=TRUE)
          spot2 <- which(!is.na(sample.spot.list[[j]]), arr.ind=TRUE)

          spot12intersect <- spot1[spot12intersect,,drop=FALSE]

          core1.center <- colMeans(apply(which(!is.na(sample.spot.core.list[[i]]), arr.ind=TRUE), 2, range))
          core2.center <- colMeans(apply(which(!is.na(sample.spot.core.list[[j]]), arr.ind=TRUE), 2, range))

          spot12assoc <- apply(spot12intersect, 1, function(x)
          {
            which.min(c(sum((core1.center - x) ^ 2), sum((core2.center - x) ^ 2)))
          })

          sample.spot.list[[j]][spot12intersect[which(spot12assoc==1),,drop=FALSE]] <- NA
          sample.spot.list[[i]][spot12intersect[which(spot12assoc==2),,drop=FALSE]] <- NA
        }
      }
    }
  }

  ## no spot detected? ##
  if( length(sample.spot.list)==0 )
  {
    util.warn("No underexpression spot detectable.")	
    sample.spot.list = list( matrix(NA, env$preferences$dim.1stLvlSom, env$preferences$dim.1stLvlSom) )
    sample.spot.list[[1]][ which.min(rowMeans(env$metadata)) ] = 1
  }
  
  ## define underexpression spots ##
  env$spot.list.underexpression <- list()

  env$spot.list.underexpression$overview.map <-
    apply(apply(env$metadata, 2, function(x) { (x - min(x)) / (max(x) - min(x)) }), 1, min)

  env$spot.list.underexpression$overview.mask <- rep(NA, env$preferences$dim.1stLvlSom ^ 2)
  env$spot.list.underexpression$spots <- list()

  for (i in seq_along(sample.spot.list))
  {
    spot.metagenes <- which(!is.na(sample.spot.list[[i]]))
    spot.genes <- rownames(env$indata)[which(env$som.result$feature.BMU %in% spot.metagenes)]

    if (length(spot.genes) > 0)
    {
      env$spot.list.underexpression$overview.mask[which(!is.na(sample.spot.list[[i]]))] <- i
      env$spot.list.underexpression$spots[[env$letters[i]]] <- list()
      env$spot.list.underexpression$spots[[env$letters[i]]]$metagenes <- spot.metagenes
      env$spot.list.underexpression$spots[[env$letters[i]]]$genes <- spot.genes
      env$spot.list.underexpression$spots[[env$letters[i]]]$mask <- rep(NA, env$preferences$dim.1stLvlSom * env$preferences$dim.1stLvlSom)
      env$spot.list.underexpression$spots[[env$letters[i]]]$mask[spot.metagenes] <- 1

      env$spot.list.underexpression$spots[[env$letters[i]]]$position <-
        colMeans(apply(env$som.result$node.summary[spot.metagenes, 1:2]+1, 2, range))

      env$spot.list.underexpression$spots[[env$letters[i]]]$beta.statistic <-
        get.beta.statistic(set.data=env$metadata[env$spot.list.underexpression$spots[[env$letters[i]]]$metagenes,,drop=FALSE],
                           weights=env$som.result$node.summary[env$spot.list.underexpression$spots[[env$letters[i]]]$metagenes,]$n.features)
    }
  }
  
  o <- order(sapply(env$spot.list.underexpression$spots, function(x)
  {
    which.min(apply(env$metadata[x$metagenes,,drop=FALSE], 2, mean))
  }))

  env$spot.list.underexpression$spots <- env$spot.list.underexpression$spots[o]
  names(env$spot.list.underexpression$spots) <- env$letters[seq_along(env$spot.list.underexpression$spots)]

  env$spot.list.underexpression$overview.mask[!is.na(env$spot.list.underexpression$overview.mask)] <-
    match(env$spot.list.underexpression$overview.mask[!is.na(env$spot.list.underexpression$overview.mask)], o)

  env$spot.list.underexpression$spotdata <-
    t(sapply(env$spot.list.underexpression$spots, function(x)
    {
      if (length(x$genes > 0))
      {
        colMeans(env$indata[x$genes,,drop=FALSE],na.rm=TRUE)
      } else
      {
        rep(0, ncol(env$indata))
      }
    }))

  colnames(env$spot.list.underexpression$spotdata) <- colnames(env$indata)

  return(env)   
}


pipeline.detectCorrelationModules <- function(env)
{
  env$spot.list.correlation <- list()
  env$spot.list.correlation$overview.map <- NA
  env$spot.list.correlation$overview.mask <- rep(NA, env$preferences$dim.1stLvlSom ^ 2)
  env$spot.list.correlation$spots <- list()

  c.map <- cor(t(env$metadata))
  diag(c.map) <- NA
  rownames(c.map) <- c(1:(env$preferences$dim.1stLvlSom * env$preferences$dim.1stLvlSom))
  colnames(c.map) <- c(1:(env$preferences$dim.1stLvlSom * env$preferences$dim.1stLvlSom))
   
  c.cluster <- rep(NA, env$preferences$dim.1stLvlSom * env$preferences$dim.1stLvlSom)
  names(c.cluster) <- c(1:(env$preferences$dim.1stLvlSom * env$preferences$dim.1stLvlSom))

  count.cluster <- 1
  
  while (is.matrix(c.map) && nrow(c.map) > 0 && ncol(c.map) > 0 && count.cluster <= 30)
  {
    start.node <- rownames(c.map)[which(c.map == max(c.map, na.rm=TRUE), arr.ind=TRUE)[1,1]]

    cluster <- names(which(c.map[start.node,] > 0.90))
    cluster <- c(start.node, cluster)

    if (length(cluster) >= env$preferences$dim.1stLvlSom / 2)
    {
      c.cluster[cluster] <- count.cluster
      geneset.genes <- rownames(env$indata)[which(env$som.result$feature.BMU %in% as.numeric(cluster))]
      
      if (length(geneset.genes) > 0)
      {
        env$spot.list.correlation$overview.mask[as.numeric(cluster)] <- count.cluster
        env$spot.list.correlation$spots[[env$LETTERS[count.cluster]]] <- list()
        env$spot.list.correlation$spots[[env$LETTERS[count.cluster]]]$metagenes <- as.numeric(cluster)
        env$spot.list.correlation$spots[[env$LETTERS[count.cluster]]]$genes <- geneset.genes
        env$spot.list.correlation$spots[[env$LETTERS[count.cluster]]]$mask <- rep(NA, env$preferences$dim.1stLvlSom * env$preferences$dim.1stLvlSom)
        env$spot.list.correlation$spots[[env$LETTERS[count.cluster]]]$mask[as.numeric(cluster)] <- 1

        env$spot.list.correlation$spots[[env$LETTERS[count.cluster]]]$position <-
          apply(apply(env$som.result$node.summary[cluster, 1:2], 2, range), 2, mean) + 0.5

        env$spot.list.correlation$spots[[env$LETTERS[count.cluster]]]$beta.statistic <-
          get.beta.statistic(set.data=env$metadata[env$spot.list.correlation$spots[[count.cluster]]$metagenes,,drop=FALSE],
                             weights=env$som.result$node.summary[env$spot.list.correlation$spots[[count.cluster]]$metagenes,]$n.features)

        count.cluster <- count.cluster + 1
      }
    }

    c.map <- c.map[-which(rownames(c.map) %in% cluster), -which(colnames(c.map) %in% cluster)]
  }

  env$spot.list.correlation$overview.map <- c.cluster
  
  o <- order(sapply(env$spot.list.correlation$spots, function(x)
  {
    which.max(apply(env$metadata[x$metagenes,,drop=FALSE], 2, mean))
  }))

  env$spot.list.correlation$spots <- env$spot.list.correlation$spots[o]
  names(env$spot.list.correlation$spots) <- env$LETTERS[seq_along(env$spot.list.correlation$spots)]

  env$spot.list.correlation$overview.mask[!is.na(env$spot.list.correlation$overview.mask)] <-
    match(env$spot.list.correlation$overview.mask[!is.na(env$spot.list.correlation$overview.mask)], sort(unique(env$spot.list.correlation$overview.mask))[o])
  env$spot.list.correlation$overview.map <- env$spot.list.correlation$overview.mask
  
  env$spot.list.correlation$spotdata <-
    t(sapply(env$spot.list.correlation$spots, function(x)
    {
      if (length(x$genes > 0))
      {
        colMeans(env$indata[x$genes,,drop=FALSE],na.rm=TRUE)
      } else
      {
        rep(0, ncol(env$indata))
      }
    }))

  colnames(env$spot.list.correlation$spotdata) <- colnames(env$indata)

  return(env)   
}


pipeline.detectKMeansModules <- function(env)
{
  n.cluster <- ceiling( env$preferences$dim.1stLvlSom / 2 )
  prototypes <- env$metadata[round(seq(1, env$preferences$dim.1stLvlSom^2, length.out=n.cluster)),]
  res <- kmeans(env$metadata, prototypes)

  env$spot.list.kmeans <- list()
  env$spot.list.kmeans$overview.map <- matrix(res$cluster, env$preferences$dim.1stLvlSom, env$preferences$dim.1stLvlSom)
  env$spot.list.kmeans$overview.mask <- rep(NA, env$preferences$dim.1stLvlSom ^ 2)
  env$spot.list.kmeans$spots <- list()

  for (i in 1:n.cluster)
  {
    nodes <- which(res$cluster == i)
    geneset.genes <- rownames(env$indata)[which(env$som.result$feature.BMU %in% nodes)]

    if (length(geneset.genes) > 0)
    {
      env$spot.list.kmeans$overview.mask[as.numeric(nodes)] <- i
      env$spot.list.kmeans$spots[[env$LETTERS[i]]] <- list()
      env$spot.list.kmeans$spots[[env$LETTERS[i]]]$metagenes <- as.numeric(nodes)
      env$spot.list.kmeans$spots[[env$LETTERS[i]]]$genes <- geneset.genes
      env$spot.list.kmeans$spots[[env$LETTERS[i]]]$mask <- rep(NA, env$preferences$dim.1stLvlSom * env$preferences$dim.1stLvlSom)
      env$spot.list.kmeans$spots[[env$LETTERS[i]]]$mask[as.numeric(nodes)] <- 1

      env$spot.list.kmeans$spots[[env$LETTERS[i]]]$position <-
        apply(apply(env$som.result$node.summary[nodes, 1:2], 2, range), 2, mean) + 0.5

      env$spot.list.kmeans$spots[[env$LETTERS[i]]]$beta.statistic <-
        get.beta.statistic(set.data=env$metadata[env$spot.list.kmeans$spots[[env$LETTERS[i]]]$metagenes,,drop=FALSE],
                           weights=env$som.result$node.summary[env$spot.list.kmeans$spots[[env$LETTERS[i]]]$metagenes,]$n.features)
    }
  }

  o <- order(sapply(env$spot.list.kmeans$spots, function(x)
  {
    which.max(apply(env$metadata[x$metagenes,,drop=FALSE], 2, mean))
  }))

  env$spot.list.kmeans$spots <- env$spot.list.kmeans$spots[o]
  names(env$spot.list.kmeans$spots) <- env$LETTERS[seq_along(env$spot.list.kmeans$spots)]
  
  env$spot.list.kmeans$overview.mask[!is.na(env$spot.list.kmeans$overview.mask)] <-
    match(env$spot.list.kmeans$overview.mask[!is.na(env$spot.list.kmeans$overview.mask)], sort(unique(env$spot.list.kmeans$overview.mask))[o])
  env$spot.list.kmeans$overview.map <- env$spot.list.kmeans$overview.mask

  env$spot.list.kmeans$spotdata <-
    t(sapply(env$spot.list.kmeans$spots, function(x)
    {
      if (length(x$genes > 0))
      {
        colMeans(env$indata[x$genes,,drop=FALSE],na.rm=TRUE)
      } else
      {
        rep(0, ncol(env$indata))
      }
    }))

  colnames(env$spot.list.kmeans$spotdata) <- colnames(env$indata)
  
  return(env)   
}


pipeline.detectGroupOverexpressionModules <- function(env)
{
  if (length(unique(env$group.labels)) > 1)
  {
    group.metadata <- do.call(cbind, by(t(env$metadata), env$group.labels, colMeans))[,unique(env$group.labels)]
    group.metadata.scaled <- apply(group.metadata, 2, function(x)   (x-min(x)) / (max(x)-min(x))   )

    ## extract sample modules ##
    sample.spot.list <- list()
    sample.spot.core.list <- list()
    n.sample.modules <- 0

    for (m in 1:ncol(group.metadata))
    {
      # define bigger core regions
      core <- matrix(NA, env$preferences$dim.1stLvlSom, env$preferences$dim.1stLvlSom)
      core[which(group.metadata.scaled[,m] > env$preferences$spot.threshold.groupmap)] <- -1

      spot.i <- 0

      while (nrow(which(core == -1, arr.ind=TRUE)) > 0)
      {
        start.pix <- which(core == -1, arr.ind=TRUE)[1,]

        spot.i <- spot.i + 1
        core <- col.pix(core, start.pix[1], start.pix[2], spot.i, env$preferences$dim.1stLvlSom)
      }
      
      # shrink each separate region to core size
      for (s.i in 1:max(core,na.rm=TRUE))
      {
        if (sum(core == s.i, na.rm=TRUE) > env$preferences$spot.coresize.groupmap)
        {
          core.metagenes <- which(core == s.i)
          o <- order(group.metadata[core.metagenes,m], decreasing=TRUE)[1:env$preferences$spot.coresize.groupmap]
          core[setdiff(core.metagenes, core.metagenes[o])] <- NA
        }
      }

      core[which(!is.na(core))] <- -1
      spot.i <- 0

      while (nrow(which(core == -1, arr.ind=TRUE)) > 0)
      {
        start.pix <- which(core == -1, arr.ind=TRUE)[1,]
        spot.i <- spot.i + 1
        core <- col.pix(core, start.pix[1], start.pix[2], spot.i, env$preferences$dim.1stLvlSom)
      }
      
      # define spot area around cores
      for (s.i in 1:max(core,na.rm=TRUE))
      {
        n.sample.modules <- n.sample.modules + 1

        sample.spot.core.list[[n.sample.modules]] <- matrix(NA, env$preferences$dim.1stLvlSom, env$preferences$dim.1stLvlSom)
        sample.spot.core.list[[n.sample.modules]][which(core == s.i)] <- env$metadata[which(core == s.i),m]

        spot <- matrix(NA, env$preferences$dim.1stLvlSom, env$preferences$dim.1stLvlSom)
        spot[which(group.metadata.scaled[,m] > env$preferences$spot.threshold.groupmap)] <- -1

        start.pix <- which(!is.na(sample.spot.core.list[[n.sample.modules]]), arr.ind=TRUE)
        start.pix <- start.pix[which.max(sample.spot.core.list[[n.sample.modules]][start.pix]),]

        spot <- col.pix(spot, start.pix[1], start.pix[2], 1, env$preferences$dim.1stLvlSom)

        sample.spot.list[[n.sample.modules]] <- matrix(NA, env$preferences$dim.1stLvlSom, env$preferences$dim.1stLvlSom)
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

    o <- order(sapply(sample.spot.core.list,function(x) mean(x,na.rm=TRUE)), decreasing=TRUE)
    sample.spot.list <- sample.spot.list[o]
    sample.spot.core.list <- sample.spot.core.list[o]

    ## merge overlapping sample cores ##
    merged <- TRUE

    while (merged)
    {
      merged <- FALSE
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

            merged <- TRUE
          } else
          {
            j <- j + 1
          }
        }

        i <- i + 1
      }
    }
     
    o <- order(sapply(sample.spot.core.list, function(x) { mean(x, na.rm=TRUE) }), decreasing=TRUE)
    sample.spot.list <- sample.spot.list[o]
    sample.spot.core.list <- sample.spot.core.list[o]

    ## shrinking overlapping spots ##
    if (length(sample.spot.list) > 1)
    {
      for (i in seq_len(length(sample.spot.list)-1))
      {
        for (j in seq(i+1, length(sample.spot.list)))
        {
          spot1 <- which(!is.na(sample.spot.list[[i]]))
          spot2 <- which(!is.na(sample.spot.list[[j]]))

          if (any(spot1 %in% spot2))
          {
            spot12intersect <- which(spot1 %in% spot2)

            spot1 <- which(!is.na(sample.spot.list[[i]]), arr.ind=TRUE)
            spot2 <- which(!is.na(sample.spot.list[[j]]), arr.ind=TRUE)

            spot12intersect <- spot1[spot12intersect, ,drop=FALSE]

            core1.center <- colMeans(apply(which(!is.na(sample.spot.core.list[[i]]), arr.ind=TRUE), 2, range))
            core2.center <- colMeans(apply(which(!is.na(sample.spot.core.list[[j]]), arr.ind=TRUE), 2, range))

            spot12assoc <- apply(spot12intersect, 1, function(x)
            {
              which.min(c(sum((core1.center - x) ^ 2), sum((core2.center - x) ^ 2)))
            })

            sample.spot.list[[j]][spot12intersect[which(spot12assoc==1),,drop=FALSE]] <- NA
            sample.spot.list[[i]][spot12intersect[which(spot12assoc==2),,drop=FALSE]] <- NA
          }
        }
      }
    }

    ## no spot detected? ##
    if( length(sample.spot.list)==0 )
    {
      util.warn("No group-overexpression spot detectable.")
      sample.spot.list = list( matrix(NA, env$preferences$dim.1stLvlSom, env$preferences$dim.1stLvlSom) )
      sample.spot.list[[1]][ which.max(rowMeans(env$metadata)) ] = 1
    }
    
    ## define overexpression spots ##
    env$spot.list.group.overexpression <- list()

    env$spot.list.group.overexpression$overview.map <-
      apply(apply(group.metadata, 2, function(x){ (x - min(x)) / (max(x) - min(x)) }), 1, max)

    env$spot.list.group.overexpression$overview.mask <- rep(NA, env$preferences$dim.1stLvlSom ^ 2)
    env$spot.list.group.overexpression$spots <- list()
     
    for (i in seq_along(sample.spot.list))
    {
      spot.metagenes <- which(!is.na(sample.spot.list[[i]]))
      spot.genes <- rownames(env$indata)[which(env$som.result$feature.BMU %in% spot.metagenes)]

      if (length(spot.genes) > 0)
      {
        env$spot.list.group.overexpression$overview.mask[which(!is.na(sample.spot.list[[i]]))] <- i
        env$spot.list.group.overexpression$spots[[env$LETTERS[i]]] <- list()
        env$spot.list.group.overexpression$spots[[env$LETTERS[i]]]$metagenes <- spot.metagenes
        env$spot.list.group.overexpression$spots[[env$LETTERS[i]]]$genes <- spot.genes
        env$spot.list.group.overexpression$spots[[env$LETTERS[i]]]$mask <- rep(NA, env$preferences$dim.1stLvlSom * env$preferences$dim.1stLvlSom)
        env$spot.list.group.overexpression$spots[[env$LETTERS[i]]]$mask[spot.metagenes] <- 1

        env$spot.list.group.overexpression$spots[[env$LETTERS[i]]]$position <-
          colMeans(apply(env$som.result$node.summary[spot.metagenes, 1:2]+1, 2, range))

        env$spot.list.group.overexpression$spots[[env$LETTERS[i]]]$beta.statistic <-
          get.beta.statistic(set.data=env$metadata[env$spot.list.group.overexpression$spots[[env$LETTERS[i]]]$metagenes,,drop=FALSE],
                             weights=env$som.result$node.summary[env$spot.list.group.overexpression$spots[[env$LETTERS[i]]]$metagenes,]$n.features)
      }
    }

    start.spot <- which.min(sapply(env$spot.list.group.overexpression$spots, function(x)
    {
      which.max(apply(group.metadata[x$metagenes,,drop=FALSE], 2, mean))
    }))

    spot.arcs <- sapply(env$spot.list.group.overexpression$spots, function(x)
    {
      -atan2(x$position['y'] - env$preferences$dim.1stLvlSom / 2, x$position['x'] - env$preferences$dim.1stLvlSom / 2)
    })

    spot.arcs <- spot.arcs - spot.arcs[start.spot]

    if (any(spot.arcs<0))
    {
      spot.arcs[which(spot.arcs<0)] <- spot.arcs[which(spot.arcs<0)] + (2 * pi)
    }

    o <- order(spot.arcs)
    
    env$spot.list.group.overexpression$spots <- env$spot.list.group.overexpression$spots[o]
    names(env$spot.list.group.overexpression$spots) <- env$LETTERS[seq_along(env$spot.list.group.overexpression$spots)]

    env$spot.list.group.overexpression$overview.mask[!is.na(env$spot.list.group.overexpression$overview.mask)] <-
      match(env$spot.list.group.overexpression$overview.mask[!is.na(env$spot.list.group.overexpression$overview.mask)], o)

    env$spot.list.group.overexpression$spotdata <-
      t(sapply(env$spot.list.group.overexpression$spots, function(x)
      {
        if (length(x$genes > 0))
        {
          colMeans(env$indata[x$genes,,drop=FALSE],na.rm=TRUE)
        } else
        {
          rep(0, ncol(env$indata))
        }
      }))

    colnames(env$spot.list.group.overexpression$spotdata) <- colnames(env$indata)
  }
  
  return(env)   
}

  
pipeline.detectDMapModules <- function(env)
{
  uh <- rep(NA, env$preferences$dim.1stLvlSom^2)
  
  for (i in 1:env$preferences$dim.1stLvlSom^2)
  {
    pos.x <- env$som.result$node.summary[i,1]
    pos.y <- env$som.result$node.summary[i,2]
    
    uh[i] <- mean(sapply(get.neighbors(pos.x, pos.y, env$preferences$dim.1stLvlSom), function(x)
    {
      x <- x - 1
      ij <- (x[1] + 1) + x[2] * env$preferences$dim.1stLvlSom
      sqrt(sum((env$metadata[ij,] - env$metadata[i,])^2))
    }))
  }
  
  uh <- matrix(log10(uh), env$preferences$dim.1stLvlSom)
  peak.matrix <- matrix( FALSE, env$preferences$dim.1stLvlSom, env$preferences$dim.1stLvlSom )

  for( pos.x in 1:env$preferences$dim.1stLvlSom )
    for( pos.y in 1:env$preferences$dim.1stLvlSom )
    {
      neighbor.dists <- sapply( get.neighbors(pos.x, pos.y, env$preferences$dim.1stLvlSom), function(x) uh[x[1],x[2]] )
      peak.matrix[pos.x, pos.y] <- all( uh[ pos.x, pos.y] <= neighbor.dists )
    }
  
  spot.matrix <- matrix(0, env$preferences$dim.1stLvlSom, env$preferences$dim.1stLvlSom )
  for( sel.minimum in which(peak.matrix)[order( uh[ which(peak.matrix) ], decreasing=TRUE )] )
  {
    spot.members <- c( sel.minimum )
    spot.d <- c( mean( uh[spot.members] ) ) 
    
    while(TRUE)
    {
      spot.neighbors <- unique( unlist(   sapply( spot.members, get.neighbors, dim=env$preferences$dim.1stLvlSom )  ) )
      spot.neighbors <- setdiff( spot.neighbors, spot.members )
      
      expansion.d <- sapply( spot.neighbors, function(x) dist( env$metadata[c(x,sel.minimum),] ) )
      expansion <- spot.neighbors[ which.min( expansion.d ) ]
      
      if( ( mean( uh[c(spot.members,expansion)] ) < spot.d[length(spot.d)] && 
            spot.d[length(spot.d)] < spot.d[length(spot.d)-1] && 
            spot.d[length(spot.d)-1] < spot.d[length(spot.d)-2] ) || 
          length(spot.members) > env$preferences$dim.1stLvlSom^2*0.1 ) 
      {
        break
      } else
      {
        spot.members <- c( spot.members, expansion )
        spot.d <- c( spot.d, mean( uh[spot.members] ) ) 
      }
      
    }
#    if( mean(uh[spot.members]) > mean(uh) ) 
    spot.matrix[spot.members] <- max(spot.matrix,na.rm=TRUE)+1
  }
  spot.matrix[which(spot.matrix==0)] <- NA
  
  
  
  env$spot.list.dmap <- list()
  env$spot.list.dmap$overview.map <- uh
  env$spot.list.dmap$overview.mask <- rep(NA, env$preferences$dim.1stLvlSom ^ 2)
  env$spot.list.dmap$spots <- list()
  
  count.cluster <- 1
  for (i in seq_along(sort(unique(na.omit(as.vector(spot.matrix))))) )
  {
    spot.metagenes <- which(spot.matrix==i)
    spot.genes <- rownames(env$indata)[which(env$som.result$feature.BMU %in% spot.metagenes)]
    
    if (length(spot.genes) > 0)
    {
      env$spot.list.dmap$overview.mask[spot.metagenes] <- count.cluster
      env$spot.list.dmap$spots[[env$LETTERS[count.cluster]]] <- list()
      env$spot.list.dmap$spots[[env$LETTERS[count.cluster]]]$metagenes <- spot.metagenes
      env$spot.list.dmap$spots[[env$LETTERS[count.cluster]]]$genes <- spot.genes
      env$spot.list.dmap$spots[[env$LETTERS[count.cluster]]]$mask <- rep(NA, env$preferences$dim.1stLvlSom * env$preferences$dim.1stLvlSom)
      env$spot.list.dmap$spots[[env$LETTERS[count.cluster]]]$mask[spot.metagenes] <- 1
      
      env$spot.list.dmap$spots[[env$LETTERS[count.cluster]]]$position <-
        colMeans(apply(env$som.result$node.summary[spot.metagenes, 1:2]+1, 2, range))
      
      env$spot.list.dmap$spots[[env$LETTERS[count.cluster]]]$beta.statistic <-
        get.beta.statistic(set.data=env$metadata[env$spot.list.dmap$spots[[env$LETTERS[count.cluster]]]$metagenes,,drop=FALSE],
                           weights=env$som.result$node.summary[env$spot.list.dmap$spots[[env$LETTERS[count.cluster]]]$metagenes,]$n.features)
      
      count.cluster <- count.cluster + 1
    }
  }
  
  env$spot.list.dmap$spotdata <-
    t(sapply(env$spot.list.dmap$spots, function(x)
    {
      if (length(x$genes > 0))
      {
        colMeans(env$indata[x$genes,,drop=FALSE],na.rm=TRUE)
      } else
      {
        rep(0, ncol(env$indata))
      }
    }))
  
  # sig.spots <- which( apply( env$spot.list.dmap$spotdata, 1, function(x) sd(x) > sd(env$spot.list.dmap$spotdata,na.rm=T) ) )
  # if( length(sig.spots) > 0 )
  # {
  #   env$spot.list.dmap$spots <- env$spot.list.dmap$spots[sig.spots]
  #   env$spot.list.dmap$overview.mask[which(!env$spot.list.dmap$overview.mask%in%sig.spots)] <- NA
  #   env$spot.list.dmap$overview.mask[!is.na(env$spot.list.dmap$overview.mask  )] <-
  #     match(env$spot.list.dmap$overview.mask[!is.na(env$spot.list.dmap$overview.mask)], sort(unique(na.omit(as.vector(env$spot.list.dmap$overview.mask)))))
  # }
  
  # start.spot <- which.min( apply( sapply(env$spot.list.dmap$spots, function(x) x$position ), 2, min ) )
  
  spot.arcs <- sapply(env$spot.list.dmap$spots, function(x)
  {
    -atan2(x$position['y'] - env$preferences$dim.1stLvlSom / 2, x$position['x'] - env$preferences$dim.1stLvlSom / 2)
  })
  
  # spot.arcs <- spot.arcs - spot.arcs[start.spot]
  
  # if (any(spot.arcs<0))
  # {
  #   spot.arcs[which(spot.arcs<0)] <- spot.arcs[which(spot.arcs<0)] + (2 * pi)
  # }
  
  
  
  o <- order(spot.arcs)
  
  env$spot.list.dmap$spots <- env$spot.list.dmap$spots[o]
  names(env$spot.list.dmap$spots) <- env$LETTERS[seq_along(env$spot.list.dmap$spots)]
  
  env$spot.list.dmap$overview.mask[!is.na(env$spot.list.dmap$overview.mask  )] <-
    match(env$spot.list.dmap$overview.mask[!is.na(env$spot.list.dmap$overview.mask)], sort(unique(na.omit(as.vector(env$spot.list.dmap$overview.mask))))[o])
  
  env$spot.list.dmap$spotdata <-
    t(sapply(env$spot.list.dmap$spots, function(x)
    {
      if (length(x$genes > 0))
      {
        colMeans(env$indata[x$genes,,drop=FALSE],na.rm=TRUE)
      } else
      {
        rep(0, ncol(env$indata))
      }
    }))
  
  colnames(env$spot.list.dmap$spotdata) <- colnames(env$indata)
  
  return(env)   
}
