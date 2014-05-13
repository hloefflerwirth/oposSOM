pipeline.genesetStatisticIntegral <- function()
{
  ### init parallel computing ###
  cl <- makeCluster(preferences$max.parallel.cores)

  ### perform GS analysis ###
  for (i in 1:length(GS.infos.overexpression$spots))
  {
    geneset.ids <- unique(na.omit(gene.ids[GS.infos.overexpression$spots[[i]]$genes]))

    GS.infos.overexpression$spots[[i]]$Fisher.p <<-
      GeneSet.Fisher(geneset.ids, unique.protein.ids, gs.def.list, sort=T, cluster=cl)
  }

  util.progress(96, 100)

  for (i in 1:length(GS.infos.underexpression$spots))
  {
    geneset.ids <- unique(na.omit(gene.ids[GS.infos.underexpression$spots[[i]]$genes]))

    GS.infos.underexpression$spots[[i]]$Fisher.p <<-
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

  for (i in 1:length(GS.infos.kmeans$spots))
  {
    geneset.ids <- unique(na.omit(gene.ids[GS.infos.kmeans$spots[[i]]$genes]))

    GS.infos.kmeans$spots[[i]]$Fisher.p <<-
      GeneSet.Fisher(geneset.ids, unique.protein.ids, gs.def.list, sort=T, cluster=cl)
  }

  util.progress(99, 100)

  if (length(unique(group.labels)) > 1)
  {
    for (i in 1:length(GS.infos.group.overexpression$spots))
    {
      geneset.ids <- unique(na.omit(gene.ids[GS.infos.group.overexpression$spots[[i]]$genes]))

      GS.infos.group.overexpression$spots[[i]]$Fisher.p <<-
        GeneSet.Fisher(geneset.ids, unique.protein.ids, gs.def.list, sort=T, cluster=cl)
    }
  }
  util.progress(100, 100)
  cat("\n")

  ### stop parallel computing ###
  try({ stopCluster(cl) }, silent=T)
}
