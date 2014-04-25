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

  if (length(GS.infos.correlation$spots) > 0)
  {
    for (i in 1:length(GS.infos.correlation$spots))
    {
      geneset.ids <- unique(na.omit(gene.ids[GS.infos.correlation$spots[[i]]$genes]))

      GS.infos.correlation$spots[[i]]$Fisher.p <<-
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

  for (i in 1:length(GS.infos.group.overexpression$spots))
  {
    geneset.ids <- unique(na.omit(gene.ids[GS.infos.group.overexpression$spots[[i]]$genes]))

    GS.infos.group.overexpression$spots[[i]]$Fisher.p <<-
      GeneSet.Fisher(geneset.ids, unique.protein.ids, gs.def.list, sort=T, cluster=cl)
  }

  util.progress(100, 100)
  util.cat("")

  ### stop parallel computing ###
  try({ stopCluster(cl) }, silent=T)
}
