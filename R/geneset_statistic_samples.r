pipeline.genesetStatisticSamples <- function()
{

  ### perform GS analysis ###
  t.ensID.m <<- t.g.m[gene.info$ensembl.mapping[,1],]
  t.ensID.m <<- do.call(rbind, by(t.ensID.m, gene.info$ensembl.mapping[,2], colMeans))
  
  mean.t.all <- colMeans( t.ensID.m )
  sd.t.all <- apply( t.ensID.m, 2, sd )
  
  if (preferences$activated.modules$geneset.analysis.exact)
  {
    gs.null.list <- list()

    for (i in seq_along(gs.def.list))
    {
      gs.null.list[[i]] <-
        list(Genes=sample(unique(gene.info$ensembl.mapping$ensembl_gene_id), length(gs.def.list[[i]]$Genes)))
    }

    null.scores <- sapply( gs.null.list, Sample.GSZ, mean.t.all, sd.t.all )
    null.culdensity <- ecdf(abs(unlist(null.scores)))
  }

  samples.GSZ.scores <<- t( sapply( gs.def.list, Sample.GSZ, mean.t.all, sd.t.all ) )
  
  spot.list.samples <<- lapply(seq_along(spot.list.samples) , function(m)
  {
    x <- spot.list.samples[[m]]

    x$GSZ.score <- samples.GSZ.scores[,i]

    if (preferences$activated.modules$geneset.analysis.exact)
    {
       x$GSZ.p.value <- 1 - null.culdensity(abs(x$GSZ.score))
       names(x$GSZ.p.value) <- names(x$GSZ.score)
    }

    return(x)
  })
  names(spot.list.samples) <<- colnames(indata)

}
