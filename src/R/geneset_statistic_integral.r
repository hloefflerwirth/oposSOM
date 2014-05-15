pipeline.genesetStatisticIntegral <- function()
{
  ### init parallel computing ###
  cl <- makeCluster(preferences$max.parallel.cores)

  ### perform GS analysis ###
  for (i in 1:length(spot.list.overexpression$spots))
  {
    geneset.ids <- unique(na.omit(gene.ids[spot.list.overexpression$spots[[i]]$genes]))

    spot.list.overexpression$spots[[i]]$Fisher.p <<-
      GeneSet.Fisher(geneset.ids, unique.protein.ids, gs.def.list, sort=T, cluster=cl)
  }

  util.progress(96, 100)

  for (i in 1:length(spot.list.underexpression$spots))
  {
    geneset.ids <- unique(na.omit(gene.ids[spot.list.underexpression$spots[[i]]$genes]))

    spot.list.underexpression$spots[[i]]$Fisher.p <<-
      GeneSet.Fisher(geneset.ids, unique.protein.ids, gs.def.list, sort=T, cluster=cl)
  }

  util.progress(97, 100)

  if (length(spot.list.correlation$spots) > 0)
  {
    for (i in 1:length(spot.list.correlation$spots))
    {
      geneset.ids <- unique(na.omit(gene.ids[spot.list.correlation$spots[[i]]$genes]))

      spot.list.correlation$spots[[i]]$Fisher.p <<-
        GeneSet.Fisher(geneset.ids, unique.protein.ids, gs.def.list, sort=T, cluster=cl)
    }
  }

  util.progress(98, 100)

  for (i in 1:length(spot.list.kmeans$spots))
  {
    geneset.ids <- unique(na.omit(gene.ids[spot.list.kmeans$spots[[i]]$genes]))

    spot.list.kmeans$spots[[i]]$Fisher.p <<-
      GeneSet.Fisher(geneset.ids, unique.protein.ids, gs.def.list, sort=T, cluster=cl)
  }

  util.progress(99, 100)

  if (length(unique(group.labels)) > 1)
  {
    for (i in 1:length(spot.list.group.overexpression$spots))
    {
      geneset.ids <- unique(na.omit(gene.ids[spot.list.group.overexpression$spots[[i]]$genes]))

      spot.list.group.overexpression$spots[[i]]$Fisher.p <<-
        GeneSet.Fisher(geneset.ids, unique.protein.ids, gs.def.list, sort=T, cluster=cl)
    }
  }

  util.progress.terminate()

  ### stop parallel computing ###
  try({ stopCluster(cl) }, silent=T)
}
