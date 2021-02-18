pipeline.genesetStatisticSamples <- function(env)
{

  ### perform GS analysis ###
  env$indata.ensID.m <- env$indata[env$gene.info$ensembl.mapping[,1],]
  env$indata.ensID.m <- do.call(rbind, by(env$indata.ensID.m, env$gene.info$ensembl.mapping[,2], colMeans))
  
  mean.ex.all <- colMeans( env$indata.ensID.m )
  sd.ex.all <- apply( env$indata.ensID.m, 2, sd )
  
  gs.null.list <- list()
  for (i in seq_along(env$gs.def.list))
  {
    gs.null.list[[i]] <-
      list(Genes=sample(unique(env$gene.info$ensembl.mapping$ensembl_gene_id), length(env$gs.def.list[[i]]$Genes)))
  }

  null.scores <- sapply( gs.null.list, Sample.GSZ, env$indata.ensID.m, mean.ex.all, sd.ex.all )
  null.culdensity <- ecdf(abs(unlist(null.scores)))

  env$samples.GSZ.scores <- t( sapply( env$gs.def.list, Sample.GSZ, env$indata.ensID.m, mean.ex.all, sd.ex.all ) )
  
  env$spot.list.samples <- lapply(seq_along(env$spot.list.samples) , function(m)
  {
    x <- env$spot.list.samples[[m]]

    x$GSZ.score <- env$samples.GSZ.scores[,m]
    x$GSZ.p.value <- 1 - null.culdensity(abs(x$GSZ.score))
    names(x$GSZ.p.value) <- names(x$GSZ.score)

    return(x)
  })
  names(env$spot.list.samples) <- colnames(env$indata)
  
  return(env)
}
