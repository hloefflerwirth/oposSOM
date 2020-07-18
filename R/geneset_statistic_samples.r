pipeline.genesetStatisticSamples <- function()
{

  ### perform GS analysis ###
  t.ensID.m <<- t.g.m[gene.info$ensembl.mapping[,1],]
  t.ensID.m <<- do.call(rbind, by(t.ensID.m, gene.info$ensembl.mapping[,2], colMeans))
  
  if (preferences$activated.modules$geneset.analysis.exact)
  {
    gs.null.list <- list()

    for (i in seq_along(gs.def.list))
    {
      gs.null.list[[i]] <-
        list(Genes=sample(unique(gene.info$ensembl.mapping$ensembl_gene_id), length(gs.def.list[[i]]$Genes)))
    }

    null.scores <- sapply( 1:ncol(indata), function(m)
    {
      all.gene.statistic <- t.ensID.m[,m]
      spot.gene.ids <- unique(gene.info$ensembl.mapping$ensembl_gene_id)

      scores <- GeneSet.GSZ(spot.gene.ids, all.gene.statistic, gs.null.list)

      return(scores)
    })

    null.culdensity <- ecdf(abs(unlist(null.scores)))
  }

  spot.list.samples <<- lapply(seq_along(spot.list.samples) , function(m)
  {
    x <- spot.list.samples[[m]]
    all.gene.statistic <- t.ensID.m[,m]
    spot.gene.ids <- unique(gene.info$ensembl.mapping$ensembl_gene_id)

    x$GSZ.score <-
      GeneSet.GSZ(spot.gene.ids, all.gene.statistic, gs.def.list, sort=FALSE)
#      GeneSet.maxmean(all.gene.statistic, gs.def.list)

    if (preferences$activated.modules$geneset.analysis.exact)
    {
       x$GSZ.p.value <- 1 - null.culdensity(abs(x$GSZ.score))
       names(x$GSZ.p.value) <- names(x$GSZ.score)
    }

    return(x)
  })
  names(spot.list.samples) <<- colnames(indata)
  
  ### GSZ table output ###
  samples.GSZ.scores <<- do.call(cbind, lapply(spot.list.samples, function(x)
  {
    return(x$GSZ.score[names(gs.def.list)])
  }))

}
