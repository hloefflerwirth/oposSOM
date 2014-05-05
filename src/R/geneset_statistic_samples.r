pipeline.genesetStatisticSamples <- function()
{
  progress.current <- 0
  progress.max <- ncol(indata)
  util.progress(progress.current, progress.max)

  ### init parallel computing ###
  cl <- parallel::makeCluster(preferences$max.parallel.cores)

  ### perform GS analysis ###
  batch.t.g.m <<- t.g.m[which(rownames(indata) %in% names(gene.ids)),]
  batch.t.g.m <<- do.call(rbind, by(batch.t.g.m, gene.ids, colMeans))

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
      all.gene.statistic <- batch.t.g.m[, m]
      spot.gene.ids <- unique.protein.ids

      null.scores <- c(null.scores, GeneSet.GSZ(spot.gene.ids, all.gene.statistic,
                                                gs.null.list, cluster=cl))

      for (spot.i in 1:length(GS.infos.samples[[m]]$spots))
      {
        spot.genes <- GS.infos.samples[[m]]$spots[[spot.i]]$genes
        spot.gene.ids <- unique(na.omit(gene.ids[spot.genes]))
        all.gene.statistic <- batch.t.g.m[, m]

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
    all.gene.statistic <- batch.t.g.m[, m]
    spot.gene.ids <- unique.protein.ids

    GS.infos.samples[[m]]$GSZ.score <<-
      GeneSet.GSZ(spot.gene.ids, all.gene.statistic, gs.def.list, sort=F, cluster=cl)

    if (preferences$geneset.analysis.exact)
    {
      GS.infos.samples[[m]]$GSZ.p.value <<-
        1 - null.culdensity(abs(GS.infos.samples[[m]]$GSZ.score))

      names(GS.infos.samples[[m]]$GSZ.p.value) <<- names(GS.infos.samples[[m]]$GSZ.score)
    }

    for (spot.i in 1:length(GS.infos.samples[[m]]$spots))
    {
      spot.genes <- GS.infos.samples[[m]]$spots[[spot.i]]$genes
      spot.gene.ids <- unique(na.omit(gene.ids[spot.genes]))

      if (length(spot.gene.ids) > 0)
      {
        GS.infos.samples[[m]]$spots[[spot.i]]$GSZ.score <<-
          GeneSet.GSZ(spot.gene.ids, all.gene.statistic, gs.def.list, sort=F, cluster=cl)

        GS.infos.samples[[m]]$spots[[spot.i]]$Fisher.p <<-
          GeneSet.Fisher(spot.gene.ids, unique.protein.ids, gs.def.list, sort=T, cluster=cl)
      } else
      {
        GS.infos.samples[[m]]$spots[[spot.i]]$GSZ.score <<- rep(0, length(gs.def.list))
        names(GS.infos.samples[[m]]$spots[[spot.i]]$GSZ.score) <<- names(gs.def.list)
        GS.infos.samples[[m]]$spots[[spot.i]]$Fisher.p <<- rep(1, length(gs.def.list))
        names(GS.infos.samples[[m]]$spots[[spot.i]]$Fisher.p) <<- names(gs.def.list)
      }

      if (preferences$geneset.analysis.exact)
      {
        GS.infos.samples[[m]]$spots[[spot.i]]$GSZ.p.value <<-
          1 - null.culdensity(abs(GS.infos.samples[[m]]$spots[[spot.i]]$GSZ.score))

        names(GS.infos.samples[[m]]$spots[[spot.i]]$GSZ.p.value) <<-
          names(GS.infos.samples[[m]]$spots[[spot.i]]$GSZ.score)
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
