pipeline.genesetStatisticSamples <- function()
{
  progress.current <- 0
  progress.max <- ncol(indata)
  util.progress(progress.current, progress.max)

  ### init parallel computing ###
  cl <- parallel::makeCluster(preferences$max.parallel.cores)

  ### perform GS analysis ###
  t.ensID.m <<- t.g.m[which(rownames(indata) %in% names(gene.ids)),]
  t.ensID.m <<- do.call(rbind, by(t.ensID.m, gene.ids, colMeans))

  if (preferences$geneset.analysis.exact)
  {
    gs.null.list <- list()

    for (i in 1:length(gs.def.list))
    {
      gs.null.list[[i]] <-
        list(Genes=sample(unique.protein.ids, length(gs.def.list[[i]]$Genes)))
    }

    null.scores <- c()

    for (m in 1:ncol(indata))
    {
      all.gene.statistic <- t.ensID.m[, m]
      spot.gene.ids <- unique.protein.ids

      null.scores <- c(null.scores, GeneSet.GSZ(spot.gene.ids, all.gene.statistic,
                                                gs.null.list, cluster=cl))

      for (spot.i in 1:length(spot.list.samples[[m]]$spots))
      {
        spot.genes <- spot.list.samples[[m]]$spots[[spot.i]]$genes
        spot.gene.ids <- unique(na.omit(gene.ids[spot.genes]))
        all.gene.statistic <- t.ensID.m[, m]

        null.scores <- c(null.scores,  GeneSet.GSZ(spot.gene.ids, all.gene.statistic,
                                                   gs.null.list, cluster=cl))
      }

      progress.current <- progress.current + 0.5
      util.progress(progress.current, progress.max)
    }

    null.culdensity <- ecdf(abs(null.scores))
  }

  for (m in 1:ncol(indata))
  {
    all.gene.statistic <- t.ensID.m[, m]
    spot.gene.ids <- unique.protein.ids

    spot.list.samples[[m]]$GSZ.score <<-
      GeneSet.GSZ(spot.gene.ids, all.gene.statistic, gs.def.list, sort=F, cluster=cl)

    if (preferences$geneset.analysis.exact)
    {
      spot.list.samples[[m]]$GSZ.p.value <<-
        1 - null.culdensity(abs(spot.list.samples[[m]]$GSZ.score))

      names(spot.list.samples[[m]]$GSZ.p.value) <<- names(spot.list.samples[[m]]$GSZ.score)
    }

    for (spot.i in 1:length(spot.list.samples[[m]]$spots))
    {
      spot.genes <- spot.list.samples[[m]]$spots[[spot.i]]$genes
      spot.gene.ids <- unique(na.omit(gene.ids[spot.genes]))

      if (length(spot.gene.ids) > 0)
      {
        spot.list.samples[[m]]$spots[[spot.i]]$GSZ.score <<-
          GeneSet.GSZ(spot.gene.ids, all.gene.statistic, gs.def.list, sort=F, cluster=cl)

        spot.list.samples[[m]]$spots[[spot.i]]$Fisher.p <<-
          GeneSet.Fisher(spot.gene.ids, unique.protein.ids, gs.def.list, sort=T, cluster=cl)
      } else
      {
        spot.list.samples[[m]]$spots[[spot.i]]$GSZ.score <<- rep(0, length(gs.def.list))
        names(spot.list.samples[[m]]$spots[[spot.i]]$GSZ.score) <<- names(gs.def.list)
        spot.list.samples[[m]]$spots[[spot.i]]$Fisher.p <<- rep(1, length(gs.def.list))
        names(spot.list.samples[[m]]$spots[[spot.i]]$Fisher.p) <<- names(gs.def.list)
      }

      if (preferences$geneset.analysis.exact)
      {
        spot.list.samples[[m]]$spots[[spot.i]]$GSZ.p.value <<-
          1 - null.culdensity(abs(spot.list.samples[[m]]$spots[[spot.i]]$GSZ.score))

        names(spot.list.samples[[m]]$spots[[spot.i]]$GSZ.p.value) <<-
          names(spot.list.samples[[m]]$spots[[spot.i]]$GSZ.score)
      }
    }

    if (preferences$geneset.analysis.exact)
    {
      progress.current <- progress.current + 0.4
    } else
    {
      progress.current <- progress.current + 0.9
    }

    util.progress(progress.current, progress.max)
  }

  ### stop parallel computing ###
  try({ stopCluster(cl) }, silent=T)
}
