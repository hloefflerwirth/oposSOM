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
  opossom$metagene.filter.list <- NULL
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
  opossom$batch.t.g.m <- NULL
  opossom$perc.DE.m <- NULL
  opossom$fdr.g.m <- NULL
  opossom$Fdr.g.m <- NULL
  opossom$all.gene.statistic <- NULL
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

  return(opossom)
}

# Executes the oposSOM pipeline.
opossom.run <- function(opossom)
{
  opossom$preferences$system.info <- Sys.info()
  opossom$preferences$started <- format(Sys.time(), "%a %b %d %X")

  # Output some info
  util.info("Started:", opossom$preferences$started)
  util.info("Setting:", opossom$preferences$dataset.name)
  util.info("1SOM Dim:", opossom$preferences$dim.som1)
  util.info("2SOM Dim:", opossom$preferences$dim.som2)

  # Prepare the environment
  environment(pipeline.prepare) <- opossom

  if (!pipeline.prepare()) {
    return()
  }

  imagename <- paste(opossom$files.name, "pre.RData")
  util.info("Saving environment image:", imagename)
  save(opossom, file=imagename)

  # Execute the pipeline
  util.info("Processing Differential Expression")
  environment(pipeline.calcStatistics) <- opossom
  pipeline.calcStatistics()

  util.info("Detecting Spots")
  environment(pipeline.detectSpotsSamples) <- opossom
  pipeline.detectSpotsSamples()

  environment(pipeline.detectSpotsIntegral) <- opossom
  pipeline.detectSpotsIntegral()

  environment(pipeline.groupAssignment) <- opossom
  pipeline.groupAssignment()

  util.info("Plotting Sample Portraits")
  environment(pipeline.sampleExpressionPortraits) <- opossom
  pipeline.sampleExpressionPortraits()

  environment(pipeline.sampleRankMaps) <- opossom
  pipeline.sampleRankMaps()

  util.info("Processing Supporting Information")
  environment(pipeline.supportingMaps) <- opossom
  pipeline.supportingMaps()

  environment(pipeline.entropyProfiles) <- opossom
  pipeline.entropyProfiles()

  environment(pipeline.topologyProfiles) <- opossom
  pipeline.topologyProfiles()

  util.info("Processing 2nd level Metagene Analysis")
  dir.create(file.path(paste(opossom$files.name, "- Results"),
                       "2nd lvl Metagene Analysis"), showWarnings=F)

  environment(pipeline.2ndLvlSimilarityAnalysis) <- opossom
  pipeline.2ndLvlSimilarityAnalysis()

  environment(pipeline.2ndLvlCorrelationAnalysis) <- opossom
  pipeline.2ndLvlCorrelationAnalysis()

  environment(pipeline.2ndLvlComponentAnalysis) <- opossom
  pipeline.2ndLvlComponentAnalysis()

  environment(pipeline.2ndLvlSom) <- opossom
  pipeline.2ndLvlSom()


  if (opossom$preferences$geneset.analysis)
  {
    util.info("Processing Geneset Analysis")
    dir.create(paste(opossom$files.name, "- Results/Geneset Analysis"),
               showWarnings=F)

    environment(pipeline.genesetStatisticSamples) <- opossom
    pipeline.genesetStatisticSamples()

    environment(pipeline.genesetStatisticIntegral) <- opossom
    pipeline.genesetStatisticIntegral()

    environment(pipeline.genesetOverviews) <- opossom
    pipeline.genesetOverviews()

    util.info("Processing Geneset Profiles and Maps")
    environment(pipeline.genesetProfilesAndMaps) <- opossom
    pipeline.genesetProfilesAndMaps()

    util.info("Processing Cancer Hallmarks")
    environment(pipeline.cancerHallmarks) <- opossom
    pipeline.cancerHallmarks()

    util.info("Processing Chromosome Expression Reports")
    environment(pipeline.chromosomeExpressionReports) <- opossom
    pipeline.chromosomeExpressionReports()
  }

  util.info("Processing Gene Lists")
  environment(pipeline.geneLists) <- opossom
  pipeline.geneLists()

  util.info("Processing Summary Sheets (Samples)")
  environment(pipeline.summarySheetsSamples) <- opossom
  pipeline.summarySheetsSamples()

  util.info("Processing Summary Sheets (Spots)")
  environment(pipeline.summarySheetsIntegral) <- opossom
  pipeline.summarySheetsIntegral()

  util.info("Processing 3rd level Spot Analysis")
  dir.create(paste(opossom$files.name, "- Results/3rd lvl Spot Analysis"),
             showWarnings=F)

  environment(pipeline.3rdLvlChromosomalEnrichment) <- opossom
  pipeline.3rdLvlChromosomalEnrichment()

  environment(pipeline.3rdLvlSummarySheets) <- opossom
  pipeline.3rdLvlSummarySheets()

  environment(pipeline.3rdLvlNetworks) <- opossom
  pipeline.3rdLvlNetworks()

  util.info("Generating HTML Report")
  environment(pipeline.htmlSummary) <- opossom
  pipeline.htmlSummary()

  environment(pipeline.htmlSampleSummary) <- opossom
  pipeline.htmlSampleSummary()

  environment(pipeline.htmlIntegralSummary) <- opossom
  pipeline.htmlIntegralSummary()

  environment(pipeline.htmlGenesetAnalysis) <- opossom
  pipeline.htmlGenesetAnalysis()

  # Save the opossom environment
  filename <- paste(opossom$files.name, ".RData", sep="")
  util.info("Saving environment image:", filename)
  save(opossom, file=filename)

  if (file.exists(imagename) && file.exists(filename))
  {
    file.remove(imagename)
  }

  # Run additional functions. (NOTE: They are changing the environment)
  environment(pipeline.groupAnalysis) <- opossom
  pipeline.groupAnalysis()

  load(filename) # Reload opossom
  environment(pipeline.differenceAnalyses) <- opossom
  pipeline.differenceAnalyses()

  load(filename) # Reload opossom
  environment(pipeline.signatureSets) <- opossom
  pipeline.signatureSets()

  util.info("Finished:", format(Sys.time(), "%a %b %d %X"))
}
