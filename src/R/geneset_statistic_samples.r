pipeline.genesetStatisticSamples <- function()
{

  ### perform GS analysis ###
  t.ensID.m <<- t.g.m[which(rownames(indata) %in% names(gene.info$ids)),]
  t.ensID.m <<- do.call(rbind, by(t.ensID.m, gene.info$ids, colMeans))

  if (preferences$activated.modules$geneset.analysis.exact)
  {
    gs.null.list <- list()

    for (i in seq_along(gs.def.list))
    {
      gs.null.list[[i]] <-
        list(Genes=sample(unique.protein.ids, length(gs.def.list[[i]]$Genes)))
    }

    null.scores <- sapply( 1:ncol(indata), function(m)
    {
      all.gene.statistic <- t.ensID.m[,m]
      spot.gene.ids <- unique.protein.ids

      scores <- GeneSet.GSZ(spot.gene.ids, all.gene.statistic, gs.null.list)

      return(scores)
    })

    null.culdensity <- ecdf(abs(unlist(null.scores)))
  }

  spot.list.samples <<- lapply(seq_along(spot.list.samples) , function(m)
  {
    x <- spot.list.samples[[m]]
    all.gene.statistic <- t.ensID.m[,m]
    spot.gene.ids <- unique.protein.ids

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
