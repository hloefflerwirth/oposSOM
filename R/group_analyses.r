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

  local.env <- new.env()
  local.env$preferences <- env$preferences
  local.env$t.ensID.m <- env$t.ensID.m
  local.env$gene.info <- env$gene.info
  local.env$gs.def.list <- env$gs.def.list
  local.env$som.result <- env$som.result
  local.env$files.name <- env$files.name
  local.env$csv.function <- env$csv.function
  local.env$color.palette.portraits <- env$color.palette.portraits
  
  # calculate differential expression statistics

  local.env$WAD.g.m <- matrix(NA, nrow(env$indata), length(unique(env$group.labels)),
                     dimnames=list(rownames(env$indata), unique(env$group.labels)))
  local.env$t.g.m <- matrix(NA, nrow(env$indata), length(unique(env$group.labels)),
                   dimnames=list(rownames(env$indata), unique(env$group.labels)))
  local.env$p.g.m <- matrix(NA, nrow(env$indata), length(unique(env$group.labels)),
                   dimnames=list(rownames(env$indata), unique(env$group.labels)))
  local.env$fdr.g.m <- matrix(NA, nrow(env$indata), length(unique(env$group.labels)),
                     dimnames=list(rownames(env$indata), unique(env$group.labels)))
  local.env$Fdr.g.m <- matrix(NA, nrow(env$indata), length(unique(env$group.labels)),
                     dimnames=list(rownames(env$indata), unique(env$group.labels)))
  local.env$n.0.m <- rep(NA, length(unique(env$group.labels)))
    names(env$n.0.m) <- unique(env$group.labels)
  local.env$perc.DE.m <- rep(NA, length(unique(env$group.labels)))
    names(env$perc.DE.m) <- unique(env$group.labels)
  

  for (gr in seq_along(unique(env$group.labels)))
  {
    samples.indata <- which(env$group.labels==unique(env$group.labels)[gr])
    
    n <- length(samples.indata)
    local.env$t.g.m[,gr] <- sqrt(n) * apply(env$indata[,samples.indata,drop=FALSE],1,function(x)
    {
      sd.estimate = if( n>1 ) sd(x) else 1
      if( sd.estimate == 0 ) sd.estimate = 1
        
      return( mean(x) / sd.estimate )
    } )
#    t.g.m[which(is.nan(t.g.m[,gr])|is.infinite(t.g.m[,gr])),gr] <- 0
    local.env$p.g.m[,gr] <- 2 - 2*pt( abs(local.env$t.g.m[,gr]), max(n-1,1) )

    suppressWarnings({
      try.res <- try({
        fdrtool.result <- fdrtool(local.env$p.g.m[,gr], statistic="pvalue", verbose=FALSE, plot=FALSE)
#        fdrtool.result <- fdrtool(t.g.m[,gr], verbose=FALSE, plot=FALSE)
      }, silent=TRUE)
    })
    
    if (!is(try.res,"try-error"))
    {
#      p.g.m[,gr] <- fdrtool.result$pval
      local.env$fdr.g.m[,gr] <- fdrtool.result$lfdr
      local.env$Fdr.g.m[,gr] <- fdrtool.result$qval
      local.env$n.0.m[gr] <- fdrtool.result$param[1,"eta0"]
      local.env$perc.DE.m[gr] <- 1 - local.env$n.0.m[gr]
    } else
    {
#      p.g.m[,gr] <- order(apply(indata[,samples.indata,drop=FALSE],1,mean)) / nrow(indata)
      local.env$fdr.g.m[,gr] <- local.env$p.g.m[,gr]
      local.env$Fdr.g.m[,gr] <- local.env$p.g.m[,gr]
      local.env$n.0.m[gr] <- 0.5
      local.env$perc.DE.m[gr] <- 0.5
    }
    
    delta.e.g.m <- apply(env$indata[,samples.indata,drop=FALSE],1,mean)

    local.env$w.g.m <- (delta.e.g.m - min(delta.e.g.m)) / (max(delta.e.g.m) - min(delta.e.g.m))
    local.env$WAD.g.m[,gr] <- local.env$w.g.m * delta.e.g.m
  }
  
  

  # average over group members
    
  local.env$metadata <- do.call(cbind, by(t(env$metadata), env$group.labels, colMeans)[unique(env$group.labels)])

  local.env$indata <- do.call(cbind, by(t(env$indata+env$indata.gene.mean),
                              env$group.labels,
                             colMeans)[unique(env$group.labels)])

  local.env$indata.gene.mean <- rowMeans(local.env$indata)

  if (local.env$preferences$feature.centralization)
  {
    local.env$indata <- local.env$indata - local.env$indata.gene.mean
  }

  local.env$group.colors <- env$group.colors[match(colnames(local.env$indata), env$group.labels)]
  local.env$group.labels <- env$group.labels[match(colnames(local.env$indata), env$group.labels)]
  names(local.env$group.labels) <- local.env$group.labels
  names(local.env$group.colors) <- local.env$group.labels




  local.env$output.paths <- c("CSV" = paste(env$files.name, "- Results/Summary Sheets - Groups/CSV Sheets"),
                     "Summary Sheets Samples"= paste(env$files.name, "- Results/Summary Sheets - Groups/Reports"))
  
  local.env <- pipeline.detectSpotsSamples(local.env)

  if (local.env$preferences$activated.modules$geneset.analysis)
  {
    local.env <- pipeline.genesetStatisticSamples(local.env)
  }

  pipeline.geneLists(local.env)
  pipeline.summarySheetsSamples(local.env)
  pipeline.htmlGroupSummary(local.env)
}
