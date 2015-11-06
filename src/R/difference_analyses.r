pipeline.differenceAnalyses = function()
{
  
  if (length(unique(group.labels)) >= 2 && length(unique(group.labels)) <= 8 )
  {
    differences.list <- apply(combn(unique(group.labels), 2), 2, function(x)
    {
      list(which(group.labels==x[1]), which(group.labels==x[2]))
    })
    
    names(differences.list) <-
      apply(combn(unique(group.labels), 2), 2, paste, collapse=" vs ")
    
  } else
  {
    differences.list <- list()
    util.warn("Skipped pairwise group analyses: too few or many groups")
  }
  
  differences.list <- c( preferences$pairwise.comparison.list, differences.list )
  
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
  
  
  util.info("Processing Differences Analyses")
  dir.create(paste(files.name, "- Results/Summary Sheets - Differences"), showWarnings=FALSE)
  dir.create(paste(files.name, "- Results/Summary Sheets - Differences/CSV Sheets"), showWarnings=FALSE)

  WAD.g.m <<- matrix(NA, nrow(indata), length(differences.list),
                     dimnames=list(rownames(indata), names(differences.list)))

  t.g.m <<- matrix(NA, nrow(indata), length(differences.list),
                   dimnames=list(rownames(indata), names(differences.list)))

  p.g.m <<- matrix(NA, nrow(indata), length(differences.list),
                   dimnames=list(rownames(indata), names(differences.list)))

  fdr.g.m <<- matrix(NA, nrow(indata), length(differences.list),
                     dimnames=list(rownames(indata), names(differences.list)))

  Fdr.g.m <<- matrix(NA, nrow(indata), length(differences.list),
                     dimnames=list(rownames(indata), names(differences.list)))

  n.0.m <<- rep(NA, length(differences.list))
  names(n.0.m) <<- names(differences.list)

  perc.DE.m <<- rep(NA, length(differences.list))
  names(perc.DE.m) <<- names(differences.list)

  indata.d <- matrix(NA, nrow(indata), length(differences.list),
                      dimnames=list(rownames(indata), names(differences.list)))
  metadata.d <- matrix(NA, nrow(metadata), length(differences.list),
                        dimnames=list(rownames(metadata), names(differences.list)))

  for (d in seq_along(differences.list))
  {
    samples.indata <-
      list(differences.list[[d]][[1]], differences.list[[d]][[2]])

    
    indata.d[,d] <- rowMeans(indata[,samples.indata[[1]],drop=FALSE]) -
                     rowMeans(indata[,samples.indata[[2]],drop=FALSE])
    
    metadata.d[,d] <- rowMeans(metadata[,samples.indata[[1]],drop=FALSE]) -
                        rowMeans(metadata[,samples.indata[[2]],drop=FALSE])
    
    
    n <- length(samples.indata[[1]])
    m <- length(samples.indata[[2]])

    S2.x <- apply(indata[,samples.indata[[1]],drop=FALSE],1,var)
    S2.x[which(S2.x==0)] <- min(S2.x[which(S2.x!=0)])
    S2.y <- apply(indata[,samples.indata[[2]],drop=FALSE],1,var)
    S2.y[which(S2.y==0)] <- min(S2.y[which(S2.y!=0)])
    
    t.g.m[,d] <<- indata.d[,d] / sqrt( S2.x/n + S2.y/m )

    df <- ( S2.x/n + S2.y/m )^2  / ( S2.x^2 / (n^2*(n-1)) + S2.y^2 / (m^2*(m-1)) )
    
    p.g.m[,d] <<- 2 - 2*pt( abs(t.g.m[,d]), df )
    
    
    suppressWarnings({
      try.res <- try({
        fdrtool.result <- fdrtool(p.g.m[,d], statistic="pvalue", verbose=FALSE, plot=FALSE)
      }, silent=TRUE)
    })

    if (class(try.res) != "try-error")
    {
      fdr.g.m[,d] <<- fdrtool.result$lfdr
      Fdr.g.m[,d] <<- fdrtool.result$qval
      n.0.m[d] <<- fdrtool.result$param[1,"eta0"]
      perc.DE.m[d] <<- 1 - n.0.m[d]
    } else
    {
      fdr.g.m[,d] <<- p.g.m[,d]
      Fdr.g.m[,d] <<- p.g.m[,d]
      n.0.m[d] <<- 0.5
      perc.DE.m[d] <<- 0.5
    }

    delta.e.g.m <- indata.d[,d]

    w.g.m <- (delta.e.g.m - min(delta.e.g.m)) / (max(delta.e.g.m) - min(delta.e.g.m))
    WAD.g.m[,d] <<- w.g.m * delta.e.g.m
  }

  indata <<- indata.d
  colnames(indata) <<- names(differences.list)

  metadata <<- metadata.d
  colnames(metadata) <<- names(differences.list)

  group.labels <<- names(differences.list)
  names(group.labels) <<- names(differences.list)
  group.colors <<- rep("gray20",length(differences.list))
  names(group.colors) <<- names(differences.list)

  output.paths <<- c("CSV" = paste(files.name, "- Results/Summary Sheets - Differences/CSV Sheets"),
                     "Summary Sheets Samples"= paste(files.name, "- Results/Summary Sheets - Differences/Reports"))

  util.call(pipeline.detectSpotsSamples, environment())

  if (preferences$geneset.analysis)
  {
    if (ncol(t.g.m) == 1)
    {
      # crack for by command, which requires >=2 columns
      t.g.m <<- cbind(t.g.m, t.g.m)
    }

    util.call(pipeline.genesetStatisticSamples, environment())
  }

  util.call(pipeline.geneLists, environment())
  util.call(pipeline.summarySheetsSamples, environment())
}
