pipeline.groupAnalysis <- function(env)
{
  dir.create(paste(env$files.name, "- Results/Summary Sheets - Groups"), showWarnings=FALSE)
  dir.create(paste(env$files.name, "- Results/Summary Sheets - Groups/CSV Sheets"), showWarnings=FALSE)

  if (env$preferences$activated.modules$geneset.analysis)
  {
    dir.create(paste(env$files.name, "- Results/Summary Sheets - Groups/Geneset Analysis"), showWarnings=FALSE)
    pipeline.groupSpecificGenesets(env)
  }

  pipeline.summarySheetsGroups(env)

  
  
  # calculate differential expression statistics

  env$WAD.g.m <- matrix(NA, nrow(env$indata), length(unique(env$group.labels)),
                     dimnames=list(rownames(env$indata), unique(env$group.labels)))
  env$t.g.m <- matrix(NA, nrow(env$indata), length(unique(env$group.labels)),
                   dimnames=list(rownames(env$indata), unique(env$group.labels)))
  env$p.g.m <- matrix(NA, nrow(env$indata), length(unique(env$group.labels)),
                   dimnames=list(rownames(env$indata), unique(env$group.labels)))
  env$fdr.g.m <- matrix(NA, nrow(env$indata), length(unique(env$group.labels)),
                     dimnames=list(rownames(env$indata), unique(env$group.labels)))
  env$Fdr.g.m <- matrix(NA, nrow(env$indata), length(unique(env$group.labels)),
                     dimnames=list(rownames(env$indata), unique(env$group.labels)))
  env$n.0.m <- rep(NA, length(unique(env$group.labels)))
    names(env$n.0.m) <- unique(env$group.labels)
  env$perc.DE.m <- rep(NA, length(unique(env$group.labels)))
    names(env$perc.DE.m) <- unique(env$group.labels)
  

  for (gr in seq_along(unique(env$group.labels)))
  {
    samples.indata <- which(env$group.labels==unique(env$group.labels)[gr])
    
    n <- length(samples.indata)
    env$t.g.m[,gr] <- sqrt(n) * apply(env$indata[,samples.indata,drop=FALSE],1,function(x)
    {
      sd.estimate = if( n>1 ) sd(x) else 1
      if( sd.estimate == 0 ) sd.estimate = 1
        
      return( mean(x) / sd.estimate )
    } )
#    t.g.m[which(is.nan(t.g.m[,gr])|is.infinite(t.g.m[,gr])),gr] <- 0
    env$p.g.m[,gr] <- 2 - 2*pt( abs(env$t.g.m[,gr]), max(n-1,1) )

    suppressWarnings({
      try.res <- try({
        fdrtool.result <- fdrtool(env$p.g.m[,gr], statistic="pvalue", verbose=FALSE, plot=FALSE)
#        fdrtool.result <- fdrtool(t.g.m[,gr], verbose=FALSE, plot=FALSE)
      }, silent=TRUE)
    })
    
    if (!is(try.res,"try-error"))
    {
#      p.g.m[,gr] <- fdrtool.result$pval
      env$fdr.g.m[,gr] <- fdrtool.result$lfdr
      env$Fdr.g.m[,gr] <- fdrtool.result$qval
      env$n.0.m[gr] <- fdrtool.result$param[1,"eta0"]
      env$perc.DE.m[gr] <- 1 - env$n.0.m[gr]
    } else
    {
#      p.g.m[,gr] <- order(apply(indata[,samples.indata,drop=FALSE],1,mean)) / nrow(indata)
      env$fdr.g.m[,gr] <- env$p.g.m[,gr]
      env$Fdr.g.m[,gr] <- env$p.g.m[,gr]
      env$n.0.m[gr] <- 0.5
      env$perc.DE.m[gr] <- 0.5
    }
    
    delta.e.g.m <- apply(env$indata[,samples.indata,drop=FALSE],1,mean)

    env$w.g.m <- (delta.e.g.m - min(delta.e.g.m)) / (max(delta.e.g.m) - min(delta.e.g.m))
    env$WAD.g.m[,gr] <- env$w.g.m * delta.e.g.m
  }
  
  

  # average over group members
    
  env$metadata <- do.call(cbind, by(t(env$metadata), env$group.labels, colMeans)[unique(env$group.labels)])

  env$indata <- do.call(cbind, by(t(env$indata+env$indata.gene.mean),
                              env$group.labels,
                             colMeans)[unique(env$group.labels)])

  env$indata.gene.mean <- rowMeans(env$indata)

  if (env$preferences$feature.centralization)
  {
    env$indata <- env$indata - env$indata.gene.mean
  }

  env$group.colors <- env$group.colors[match(colnames(env$indata), env$group.labels)]
  env$group.labels <- env$group.labels[match(colnames(env$indata), env$group.labels)]
  names(env$group.labels) <- env$group.labels
  names(env$group.colors) <- env$group.labels




  output.paths <- c("CSV" = paste(env$files.name, "- Results/Summary Sheets - Groups/CSV Sheets"),
                     "Summary Sheets Samples"= paste(env$files.name, "- Results/Summary Sheets - Groups/Reports"))
  
  env <- pipeline.detectSpotsSamples(env)

  if (env$preferences$activated.modules$geneset.analysis)
  {
    env <- pipeline.genesetStatisticSamples(env)
  }

  pipeline.geneLists(env)
  pipeline.summarySheetsSamples(env)
  pipeline.htmlGroupSummary(env)
}
