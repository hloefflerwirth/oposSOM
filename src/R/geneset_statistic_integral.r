pipeline.genesetStatisticIntegral <- function()
{
  spot.fisher.p <- function(spot)
  {
    spot$Fisher.p <- GeneSet.Fisher(unique(na.omit(gene.ids[spot$genes])),
                                    unique.protein.ids, gs.def.list, sort=TRUE)

    return(spot)
  }

  spot.list.overexpression$spots <<- lapply( spot.list.overexpression$spots, spot.fisher.p)

  util.progress(96, 100)

  spot.list.underexpression$spots <<- lapply( spot.list.underexpression$spots, spot.fisher.p)

  util.progress(97, 100)

  if (length(spot.list.correlation$spots) > 0)
  {
    spot.list.correlation$spots <<- lapply( spot.list.correlation$spots, spot.fisher.p)
  }

  util.progress(98, 100)

  spot.list.kmeans$spots <<- lapply( spot.list.kmeans$spots, spot.fisher.p)

  util.progress(99, 100)

  if (length(unique(group.labels)) > 1)
  {
    spot.list.group.overexpression$spots <<- lapply( spot.list.group.overexpression$spots, spot.fisher.p)
  }
  
  spot.list.dmap$spots <<- lapply( spot.list.dmap$spots, spot.fisher.p)

  util.progress.terminate()

}
