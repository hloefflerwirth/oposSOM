pipeline.genesetStatisticSamples <- function(env)
{
  ### perform GS analysis ###
  env$indata.ensID.m <- env$indata[env$gene.info$ensembl.mapping[,1],]
  env$indata.ensID.m <- do.call(rbind, by(env$indata.ensID.m, env$gene.info$ensembl.mapping[,2], colMeans))
  
  mean.ex.all <- colMeans( env$indata.ensID.m )
  sd.ex.all <- apply( env$indata.ensID.m, 2, sd )
  
  env$samples.GSZ.scores <- t( sapply( env$gs.def.list, Sample.GSZ, env$indata.ensID.m, mean.ex.all, sd.ex.all ) )

  return(env)
}
