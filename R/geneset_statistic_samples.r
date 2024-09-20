pipeline.genesetStatisticSamples <- function(env)
{
  
  if( is.null(env$indata.ensID.m) )
  {
    env$indata.ensID.m <- env$indata[env$gene.info$ensembl.mapping[,1],]
    env$indata.ensID.m <- do.call(rbind, by.minicluster(env$indata.ensID.m, env$gene.info$ensembl.mapping[,2], colMeans))
  }

  mean.ex.all <- colMeans( env$indata.ensID.m )
  sd.ex.all <- apply( env$indata.ensID.m, 2, sd )

  env$samples.GSZ.scores <- t( chunk.sapply( env$gs.def.list, function(gene.set)
  {
    mean.ex.set <- colMeans( ex.ensID.m[gene.set$Genes,] )
    
    GSZ <- sqrt(length(gene.set$Genes)) * ( mean.ex.set - mean.ex.all ) / sd.ex.all
    
  }, list(ex.ensID.m=env$indata.ensID.m, mean.ex.all=mean.ex.all, sd.ex.all=sd.ex.all ), max.chunks=4 ) )
  
  
  return(env)
}
