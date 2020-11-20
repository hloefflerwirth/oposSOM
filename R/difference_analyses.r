pipeline.differenceAnalyses = function(env)
{
  
  if (length(unique(env$group.labels)) >= 2 && length(unique(env$group.labels)) <= 8 )
  {
    differences.list <- apply(combn(unique(env$group.labels), 2), 2, function(x)
    {
      list(which(env$group.labels==x[1]), which(env$group.labels==x[2]))
    })
    
    names(differences.list) <-
      apply(combn(unique(env$group.labels), 2), 2, paste, collapse=" vs ")
    
  } else
  {
    differences.list <- list()
    util.warn("Skipped pairwise group analyses: too few or many groups")
  }
  
  differences.list <- c( env$preferences$pairwise.comparison.list, differences.list )
  
  singleton.differences <- sapply( differences.list, function(x) length(x[[1]])<2 || length(x[[2]])<2 )
  if( any(singleton.differences) )
  {
    differences.list <- differences.list[which(!singleton.differences)]
    util.warn("Skipped difference analysis for groups with only one sample")
  }
  
  if (length(differences.list) == 0)
  {
    return()
  }
  
  
  dir.create(paste(env$files.name, "- Results/Summary Sheets - Differences"), showWarnings=FALSE)
  dir.create(paste(env$files.name, "- Results/Summary Sheets - Differences/CSV Sheets"), showWarnings=FALSE)

  env$WAD.g.m <- matrix(NA, nrow(env$indata), length(differences.list),
                     dimnames=list(rownames(env$indata), names(differences.list)))

  env$t.g.m <- matrix(NA, nrow(env$indata), length(differences.list),
                   dimnames=list(rownames(env$indata), names(differences.list)))

  env$p.g.m <- matrix(NA, nrow(env$indata), length(differences.list),
                   dimnames=list(rownames(env$indata), names(differences.list)))

  env$fdr.g.m <- matrix(NA, nrow(env$indata), length(differences.list),
                     dimnames=list(rownames(env$indata), names(differences.list)))

  env$Fdr.g.m <- matrix(NA, nrow(env$indata), length(differences.list),
                     dimnames=list(rownames(env$indata), names(differences.list)))

  n.0.m <- rep(NA, length(differences.list))
  names(n.0.m) <- names(differences.list)

  perc.DE.m <- rep(NA, length(differences.list))
  names(perc.DE.m) <- names(differences.list)

  indata.d <- matrix(NA, nrow(env$indata), length(differences.list),
                      dimnames=list(rownames(env$indata), names(differences.list)))
  metadata.d <- matrix(NA, nrow(env$metadata), length(differences.list),
                        dimnames=list(rownames(env$metadata), names(differences.list)))

  for (d in seq_along(differences.list))
  {
    samples.indata <-
      list(differences.list[[d]][[1]], differences.list[[d]][[2]])

    
    indata.d[,d] <- rowMeans(env$indata[,samples.indata[[1]],drop=FALSE]) -
                     rowMeans(env$indata[,samples.indata[[2]],drop=FALSE])
    
    metadata.d[,d] <- rowMeans(env$metadata[,samples.indata[[1]],drop=FALSE]) -
                        rowMeans(env$metadata[,samples.indata[[2]],drop=FALSE])
    
    
    n <- length(samples.indata[[1]])
    m <- length(samples.indata[[2]])

    S2.x <- apply(env$indata[,samples.indata[[1]],drop=FALSE],1,var)
    S2.x[which(S2.x==0)] <- min(S2.x[which(S2.x!=0)])
    S2.y <- apply(env$indata[,samples.indata[[2]],drop=FALSE],1,var)
    S2.y[which(S2.y==0)] <- min(S2.y[which(S2.y!=0)])
    
    env$t.g.m[,d] <- indata.d[,d] / sqrt( S2.x/n + S2.y/m )

    df <- ( S2.x/n + S2.y/m )^2  / ( S2.x^2 / (n^2*(n-1)) + S2.y^2 / (m^2*(m-1)) )
    
    env$p.g.m[,d] <- 2 - 2*pt( abs(env$t.g.m[,d]), df )
    
    
    suppressWarnings({
      try.res <- try({
        fdrtool.result <- fdrtool(env$p.g.m[,d], statistic="pvalue", verbose=FALSE, plot=FALSE)
      }, silent=TRUE)
    })

    if (!is(try.res,"try-error"))
    {
      env$fdr.g.m[,d] <- fdrtool.result$lfdr
      env$Fdr.g.m[,d] <- fdrtool.result$qval
      n.0.m[d] <- fdrtool.result$param[1,"eta0"]
      perc.DE.m[d] <- 1 - n.0.m[d]
    } else
    {
      env$fdr.g.m[,d] <- env$p.g.m[,d]
      env$Fdr.g.m[,d] <- env$p.g.m[,d]
      n.0.m[d] <- 0.5
      perc.DE.m[d] <- 0.5
    }

    delta.e.g.m <- indata.d[,d]

    w.g.m <- (delta.e.g.m - min(delta.e.g.m)) / (max(delta.e.g.m) - min(delta.e.g.m))
    env$WAD.g.m[,d] <- w.g.m * delta.e.g.m
  }

  indata <- indata.d
  colnames(indata) <- names(differences.list)

  metadata <- metadata.d
  colnames(metadata) <- names(differences.list)

  group.labels <- names(differences.list)
  names(group.labels) <- names(differences.list)
  group.colors <- rep("gray20",length(differences.list))
  names(group.colors) <- names(differences.list)

  output.paths <- c("CSV" = paste(env$files.name, "- Results/Summary Sheets - Differences/CSV Sheets"),
                     "Summary Sheets Samples"= paste(env$files.name, "- Results/Summary Sheets - Differences/Reports"))

  env <- pipeline.detectSpotsSamples(env)

  if (env$preferences$activated.modules$geneset.analysis)
  {
    if (ncol(env$t.g.m) == 1)
    {
      # crack for by command, which requires >=2 columns
      env$t.g.m <- cbind(env$t.g.m, env$t.g.m)
    }

    env <- pipeline.genesetStatisticSamples(env)
  }

  pipeline.geneLists(env)
  pipeline.summarySheetsSamples(env)
  pipeline.htmlDifferencesSummary(env)
}
