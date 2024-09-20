pipeline.groupAnalysis <- function(env)
{
  dir.create("Summary Sheets - Groups", showWarnings=FALSE)
  dir.create("Summary Sheets - Groups/CSV Sheets", showWarnings=FALSE)

  if (env$preferences$activated.modules$geneset.analysis)
  {
    dir.create("Summary Sheets - Groups/Geneset Analysis", showWarnings=FALSE)
    pipeline.groupSpecificGenesets(env)
  }

  pipeline.summarySheetsGroups(env)

  local.env <- new.env()
  local.env$preferences <- env$preferences
  local.env$gene.info <- env$gene.info
  local.env$gs.def.list <- env$gs.def.list
  local.env$som.result <- env$som.result
  local.env$files.name <- env$files.name
  local.env$csv.function <- env$csv.function
  local.env$color.palette.portraits <- env$color.palette.portraits
  local.env$color.palette.heatmaps <- env$color.palette.heatmaps
  
  # calculate differential expression statistics

  local.env$p.g.m <- matrix(NA, nrow(env$indata), length(unique(env$group.labels)),
                   dimnames=list(rownames(env$indata), unique(env$group.labels)))
  local.env$fdr.g.m <- matrix(NA, nrow(env$indata), length(unique(env$group.labels)),
                     dimnames=list(rownames(env$indata), unique(env$group.labels)))
  local.env$n.0.m <- rep(NA, length(unique(env$group.labels)))
  names(local.env$n.0.m) <- unique(env$group.labels)
  local.env$perc.DE.m <- rep(NA, length(unique(env$group.labels)))
  names(local.env$perc.DE.m) <- unique(env$group.labels)
  

  for (gr in seq_along(unique(env$group.labels)))
  {
    samples.indata <- which(env$group.labels==unique(env$group.labels)[gr])

    local.env$p.g.m[,gr] <- apply( env$indata, 1, function(x)
    {
      if( length(x[-samples.indata])<2 || var(x[-samples.indata]) == 0 ) return(1) 
      p <- t.test( x[samples.indata], x[-samples.indata], var.equal=length(samples.indata)==1 )$p.value
			if( p < 1e-16) p <- 1e-16
      return( p )
    } )
      
    suppressWarnings({
      try.res <- try({
        fdrtool.result <- fdrtool(local.env$p.g.m[,gr], statistic="pvalue", verbose=FALSE, plot=FALSE)
      }, silent=TRUE)
    })
    
    if (!is(try.res,"try-error"))
    {
      local.env$fdr.g.m[,gr] <- fdrtool.result$lfdr
      local.env$n.0.m[gr] <- fdrtool.result$param[1,"eta0"]
      local.env$perc.DE.m[gr] <- 1 - local.env$n.0.m[gr]
    } else
    {
      local.env$fdr.g.m[,gr] <- local.env$p.g.m[,gr]
      local.env$n.0.m[gr] <- 0.5
      local.env$perc.DE.m[gr] <- 0.5
    }    
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



  local.env$output.paths <- c("CSV" = "Summary Sheets - Groups/CSV Sheets",
                     "Summary Sheets Samples"= "Summary Sheets - Groups/Reports")
  
  if (local.env$preferences$activated.modules$geneset.analysis)
  {
    local.env <- pipeline.genesetStatisticSamples(local.env)
  }

  pipeline.summarySheetsSamples(local.env)
  pipeline.htmlGroupSummary(local.env)
}
