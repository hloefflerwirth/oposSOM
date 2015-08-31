pipeline.groupAnalysis <- function()
{
  if (length(unique(group.labels)) < 2)
  {
    return()
  }

  util.info("Processing Group-centered Analyses")
  dir.create(paste(files.name, "- Results/Summary Sheets - Groups"), showWarnings=FALSE)
  dir.create(paste(files.name, "- Results/Summary Sheets - Groups/CSV Sheets"), showWarnings=FALSE)

  if (preferences$geneset.analysis)
  {
    dir.create(paste(files.name, "- Results/Summary Sheets - Groups/Geneset Analysis"), showWarnings=FALSE)
    util.call(pipeline.groupSpecificGenesets, environment())
  }

  util.call(pipeline.summarySheetsGroups, environment())

  
  
  # calculate differential expression statistics

  WAD.g.m <<- matrix(NA, nrow(indata), length(unique(group.labels)),
                     dimnames=list(rownames(indata), unique(group.labels)))
  t.g.m <<- matrix(NA, nrow(indata), length(unique(group.labels)),
                   dimnames=list(rownames(indata), unique(group.labels)))
  p.g.m <<- matrix(NA, nrow(indata), length(unique(group.labels)),
                   dimnames=list(rownames(indata), unique(group.labels)))
  fdr.g.m <<- matrix(NA, nrow(indata), length(unique(group.labels)),
                     dimnames=list(rownames(indata), unique(group.labels)))
  Fdr.g.m <<- matrix(NA, nrow(indata), length(unique(group.labels)),
                     dimnames=list(rownames(indata), unique(group.labels)))
  n.0.m <<- rep(NA, length(unique(group.labels)))
    names(n.0.m) <<- unique(group.labels)
  perc.DE.m <<- rep(NA, length(unique(group.labels)))
    names(perc.DE.m) <<- unique(group.labels)
  

  for (gr in seq_along(unique(group.labels)))
  {
    samples.indata <- which(group.labels==unique(group.labels)[gr])
    
    n <- length(samples.indata)
    t.g.m[,gr] <<- sqrt(n) * apply(indata[,samples.indata,drop=FALSE],1,function(x) mean(x) / sd(x) )
    p.g.m[,gr] <<- 2 - 2*pt( abs(t.g.m[,gr]), n-1 )

    suppressWarnings({
      try.res <- try({
        fdrtool.result <- fdrtool(p.g.m[,gr], statistic="pvalue", verbose=FALSE, plot=FALSE)
#        fdrtool.result <- fdrtool(t.g.m[,gr], verbose=FALSE, plot=FALSE)
      }, silent=TRUE)
    })
    
    if (class(try.res) != "try-error")
    {
#      p.g.m[,gr] <<- fdrtool.result$pval
      fdr.g.m[,gr] <<- fdrtool.result$lfdr
      Fdr.g.m[,gr] <<- fdrtool.result$qval
      n.0.m[gr] <<- fdrtool.result$param[1,"eta0"]
      perc.DE.m[gr] <<- 1 - n.0.m[gr]
    } else
    {
#      p.g.m[,gr] <<- order(apply(indata[,samples.indata,drop=FALSE],1,mean)) / nrow(indata)
      fdr.g.m[,gr] <<- p.g.m[,gr]
      Fdr.g.m[,gr] <<- p.g.m[,gr]
      n.0.m[gr] <<- 0.5
      perc.DE.m[gr] <<- 0.5
    }
    
    delta.e.g.m <- apply(indata[,samples.indata,drop=FALSE],1,mean)

    w.g.m <- (delta.e.g.m - min(delta.e.g.m)) / (max(delta.e.g.m) - min(delta.e.g.m))
    WAD.g.m[,gr] <<- w.g.m * delta.e.g.m
  }
  
  

  # average over group members
    
  metadata <<- do.call(cbind, by(t(metadata), group.labels, colMeans)[unique(group.labels)])

  indata <<- do.call(cbind, by(t(indata+indata.gene.mean),
                             group.labels,
                             colMeans)[unique(group.labels)])

  indata.gene.mean <<- rowMeans(indata)

  if (preferences$feature.centralization)
  {
    indata <<- indata - indata.gene.mean
  }

  group.colors <<- group.colors[match(colnames(indata), group.labels)]
  group.labels <<- group.labels[match(colnames(indata), group.labels)]
  names(group.labels) <<- group.labels
  names(group.colors) <<- group.labels




  output.paths <<- c("CSV" = paste(files.name, "- Results/Summary Sheets - Groups/CSV Sheets"),
                     "Summary Sheets Samples"= paste(files.name, "- Results/Summary Sheets - Groups/Reports"))
  
  util.call(pipeline.detectSpotsSamples, environment())

  if (preferences$geneset.analysis)
  {
    util.call(pipeline.genesetStatisticSamples, environment())
  }

  util.call(pipeline.geneLists, environment())
  util.call(pipeline.summarySheetsSamples, environment())
}
