# Creates a new opossom environment
opossom.new <- function(preferences)
{
  # Load external packages
  require.cran("som")
  require.cran("fastICA")
  require.cran("scatterplot3d")
  require.cran("pixmap")
  require.cran("fdrtool")
  require.cran("igraph")
  require.cran("ape")
  require.cran("KernSmooth")
  require.cran("parallel")
  require.cran("foreach")

  if (Sys.info()["sysname"] == "Windows") {
    require.cran("doSNOW")
  } else {
    require.cran("doMC")
  }

  # Init the environment
  opossom <- new.env()
  opossom$indata <- NULL
  opossom$indata.original <- NULL
  opossom$indata.sample.mean <- NULL
  opossom$indata.sample.var <- NULL
  opossom$indata.mean.level <- NULL
  opossom$group.labels <- NULL
  opossom$group.colors <- NULL
  opossom$group.metadata <- NULL
  opossom$group.bootstrap.score <- NULL
  opossom$unique.group.colors <- NULL
  opossom$files.name <- NULL
  opossom$output.paths <- NULL
  opossom$colramp <- NULL
  opossom$gene.ids <- NULL
  opossom$gene.names <- NULL
  opossom$gene.descriptions <- NULL
  opossom$gene.positions <- NULL
  opossom$gene.positions.list <- NULL
  opossom$gene.positions.table <- NULL
  opossom$gs.def.list <- NULL
  opossom$gs.def.list.categories <- NULL
  opossom$unique.protein.ids <- NULL
  opossom$som.result <- NULL
  opossom$som.nodes <- NULL
  opossom$metadata <- NULL
  opossom$loglog.metadata <- NULL
  opossom$loglog.group.metadata <- NULL
  opossom$WAD.metadata <- NULL
  opossom$WAD.group.metadata <- NULL
  opossom$genes.coordinates <- NULL
  opossom$mean.e.g <- NULL
  opossom$e.g.m <- NULL
  opossom$delta.e.g.m <- NULL
  opossom$sd.g.m <- NULL
  opossom$LPE.g.m <- NULL
  opossom$WAD.g.m <- NULL
  opossom$t.g.m <- NULL
  opossom$p.g.m <- NULL
  opossom$w.g.m <- NULL
  opossom$n.0.m <- NULL
  opossom$t.m <- NULL
  opossom$p.m <- NULL
  opossom$fdr.m <- NULL
  opossom$perc.DE.m <- NULL
  opossom$fdr.g.m <- NULL
  opossom$Fdr.g.m <- NULL
  opossom$mean.LPE2 <- NULL
  opossom$GS.infos.samples <- NULL
  opossom$GS.infos.overexpression <- NULL
  opossom$GS.infos.underexpression <- NULL
  opossom$GS.infos.correlation <- NULL
  opossom$GS.infos.kmeans <- NULL
  opossom$GS.infos.group.overexpression <- NULL

  # Generate some additional letters
  opossom$LETTERS <- c(LETTERS, as.vector(sapply(1:10, function(x) {
    paste(LETTERS, x, sep="")
  })))

  opossom$letters <- c(letters, as.vector(sapply(1:10, function(x) {
    paste(letters, x, sep="")
  })))

  # Set default preferences
  opossom$preferences <- list(dataset.name = "Unnamed",
                              error.model = "all.samples",
                              dim.som1 = 20,
                              dim.som2 = 20,
                              training.extension = 1,
                              rotate.som1 = 0,
                              flip.som1 = F,
                              ensembl.dataset = "",
                              ensembl.rowname.ids = "",
                              geneset.analysis = F,
                              geneset.analysis.exact = F,
                              max.parallel.cores = detectCores() / 2,
                              sample.spot.cutoff = 0.65,
                              summary.spot.core = 3,
                              summary.spot.threshold = 0.95,
                              group.spot.core = 5,
                              group.spot.threshold = 0.75,
                              feature.mean.normalization = T,
                              sample.quantile.normalization = T,
                              differences.list = list())

  # Merge user supplied preferences
  for (key in intersect(names(opossom$preferences), names(preferences))) {
    opossom$preferences[key] <- preferences[key]
  }

  if (!opossom$preferences$error.model %in% c("replicates", "all.samples", "groups")) {
    opossom$preferences$error.model <- "all.samples"
    util.warn("Invalid value of \"error.model\". Using \"all.samples\"")
  }

  # Load opossom genesets
  data(opossom.genesets)
  opossom$preferences$geneset.custom.list <- opossom.genesets

  return(opossom)
}

# Executes the oposSOM pipeline.
opossom.run <- function(opossom)
{
  opossom$preferences$system.info <- Sys.info()
  opossom$preferences$started <- format(Sys.time(), "%a %b %d %X\n")
  environment(pipeline.run) <- opossom
  pipeline.run()
}
