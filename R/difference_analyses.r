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
  
  local.env <- new.env()
  local.env$preferences <- env$preferences
  local.env$gene.info <- env$gene.info
  local.env$gs.def.list <- env$gs.def.list
  local.env$som.result <- env$som.result
  local.env$files.name <- env$files.name
  local.env$csv.function <- env$csv.function
  local.env$color.palette.portraits <- env$color.palette.portraits
  local.env$color.palette.heatmaps <- env$color.palette.heatmaps
  local.env$indata.gene.mean <- env$indata.gene.mean

  local.env$p.g.m <- matrix(NA, nrow(env$indata), length(differences.list),
                   dimnames=list(rownames(env$indata), names(differences.list)))

  local.env$fdr.g.m <- matrix(NA, nrow(env$indata), length(differences.list),
                     dimnames=list(rownames(env$indata), names(differences.list)))

  local.env$n.0.m <- rep(NA, length(differences.list))
  names(local.env$n.0.m) <- names(differences.list)

  local.env$perc.DE.m <- rep(NA, length(differences.list))
  names(local.env$perc.DE.m) <- names(differences.list)

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
    
    local.env$p.g.m[,d] <- apply( env$indata, 1, function(x)
    {
      if( length(samples.indata[[1]])>1 && var(x[samples.indata[[1]]]) == 0 ) return(1) 
      if( length(samples.indata[[2]])>1 && var(x[samples.indata[[2]]]) == 0 ) return(1) 
      
      return( t.test( x[samples.indata[[1]]], x[samples.indata[[2]]], paired=FALSE, var.equal=(length(samples.indata[[1]])==1 || length(samples.indata[[2]])==1 ) )$p.value )
    } )
    
    suppressWarnings({
      try.res <- try({
        fdrtool.result <- fdrtool(local.env$p.g.m[,d], statistic="pvalue", verbose=FALSE, plot=FALSE)
      }, silent=TRUE)
    })

    if (!is(try.res,"try-error"))
    {
      local.env$fdr.g.m[,d] <- fdrtool.result$lfdr
      local.env$n.0.m[d] <- fdrtool.result$param[1,"eta0"]
      local.env$perc.DE.m[d] <- 1 - local.env$n.0.m[d]
    } else
    {
      local.env$fdr.g.m[,d] <- local.env$p.g.m[,d]
      local.env$n.0.m[d] <- 0.5
      local.env$perc.DE.m[d] <- 0.5
    }

  }

  local.env$indata <- indata.d
  colnames(local.env$indata) <- names(differences.list)

  local.env$metadata <- metadata.d
  colnames(local.env$metadata) <- names(differences.list)

  local.env$group.labels <- names(differences.list)
  names(local.env$group.labels) <- names(differences.list)
  local.env$group.colors <- rep("gray20",length(differences.list))
  names(local.env$group.colors) <- names(differences.list)

  local.env$output.paths <- c("CSV" = paste(env$files.name, "- Results/Summary Sheets - Differences/CSV Sheets"),
                     "Summary Sheets Samples"= paste(env$files.name, "- Results/Summary Sheets - Differences/Reports"))

  local.env <- pipeline.detectSpotsSamples(local.env)

  if (local.env$preferences$activated.modules$geneset.analysis)
  {
    if (ncol(local.env$indata) == 1)   # crack for by command, which requires >=2 columns
    {
      local.env$indata <- cbind(local.env$indata, local.env$indata)
      local.env <- pipeline.detectSpotsSamples(local.env)
      local.env <- pipeline.genesetStatisticSamples(local.env)
      local.env$indata <- local.env$indata[,1,drop=FALSE]
      local.env$spot.list.samples <- local.env$spot.list.samples[1]
      local.env$samples.GSZ.scores <- local.env$samples.GSZ.scores[,1,drop=FALSE]
    } else
    {
      local.env <- pipeline.genesetStatisticSamples(local.env)
    }
  }

  pipeline.geneLists(local.env)
  pipeline.summarySheetsSamples(local.env)
  pipeline.htmlDifferencesSummary(local.env)
}
