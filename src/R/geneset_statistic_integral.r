pipeline.genesetStatisticIntegral <- function()
{
  spot.fisher.p <- function(spot)
  {
    spot$Fisher.p <- GeneSet.Fisher(unique(na.omit(gene.ids[spot$genes])),
                                    unique.protein.ids, gs.def.list, sort=TRUE)

    return(spot)
  }

  ### init parallel computing ###
  cl <- makeCluster(preferences$max.parallel.cores)

  spot.list.overexpression$spots <<-
    parLapply(cl, spot.list.overexpression$spots, spot.fisher.p)

  util.progress(96, 100)

  spot.list.underexpression$spots <<-
    parLapply(cl, spot.list.underexpression$spots, spot.fisher.p)

  util.progress(97, 100)

  if (length(spot.list.correlation$spots) > 0)
  {
    spot.list.correlation$spots <<-
      parLapply(cl, spot.list.correlation$spots, spot.fisher.p)
  }

  util.progress(98, 100)

  spot.list.kmeans$spots <<-
    parLapply(cl, spot.list.kmeans$spots, spot.fisher.p)

  util.progress(99, 100)

  if (length(unique(group.labels)) > 1)
  {
    spot.list.group.overexpression$spots <<-
      parLapply(cl, spot.list.group.overexpression$spots, spot.fisher.p)
  }

  util.progress.terminate()

  ### stop parallel computing ###
  try({ stopCluster(cl) }, silent=TRUE)
}
