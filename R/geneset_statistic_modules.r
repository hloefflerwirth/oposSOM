pipeline.genesetStatisticModules <- function(env)
{
  plan(multisession, workers = min(8,env$preferences$max.cores))
  oopts <- options(future.globals.maxSize = 4.0 * 1e9)
  on.exit(options(oopts))
  
  
  ### perform GS enrichment analysis ###
  
  gene.info <- env$gene.info
  gs.def.list <- env$gs.def.list
  
  spot.fisher.p <- function(spot)
  {
    spot$Fisher.p <- GeneSet.Fisher(unique(gene.info$ensembl.mapping$ensembl_gene_id[ which(gene.info$ensembl.mapping[,1]%in%spot$genes) ]),
                                    unique(gene.info$ensembl.mapping$ensembl_gene_id), gs.def.list, sort=TRUE)

    return(spot)
  }

  env$spot.list.overexpression$spots <- future_lapply( env$spot.list.overexpression$spots, spot.fisher.p)
  env$spot.list.underexpression$spots <- future_lapply( env$spot.list.underexpression$spots, spot.fisher.p)
  if (length(env$spot.list.correlation$spots) > 0)
  {
    env$spot.list.correlation$spots <- future_lapply( env$spot.list.correlation$spots, spot.fisher.p)
  }
  env$spot.list.kmeans$spots <- future_lapply( env$spot.list.kmeans$spots, spot.fisher.p)
  if (length(unique(env$group.labels)) > 1)
  {
    env$spot.list.group.overexpression$spots <- future_lapply( env$spot.list.group.overexpression$spots, spot.fisher.p)
  }  
  env$spot.list.dmap$spots <- future_lapply( env$spot.list.dmap$spots, spot.fisher.p)

  return(env)
}
