# Creates a new opossom environment
opossom.new <- function(preferences)
{
  # Init the environment
  opossom <- new.env()
  opossom$t.ensID.m <- NULL
  opossom$colramp <- NULL
  opossom$Fdr.g.m <- NULL
  opossom$fdr.g.m <- NULL
  opossom$fdr.m <- NULL
  opossom$files.name <- NULL
  opossom$gene.descriptions <- NULL
  opossom$gene.ids <- NULL
  opossom$gene.names <- NULL
  opossom$gene.positions <- NULL
  opossom$gene.positions.list <- NULL
  opossom$gene.positions.table <- NULL
  opossom$gene.coordinates <- NULL
  opossom$group.bootstrap.score <- NULL
  opossom$group.colors <- NULL
  opossom$group.labels <- NULL
  opossom$group.metadata <- NULL
  opossom$gs.def.list <- NULL
  opossom$gs.def.list.categories <- NULL
  opossom$spot.list.correlation <- NULL
  opossom$spot.list.group.overexpression <- NULL
  opossom$spot.list.kmeans <- NULL
  opossom$spot.list.overexpression <- NULL
  opossom$spot.list.samples <- NULL
  opossom$spot.list.underexpression <- NULL
  opossom$indata <- NULL
  opossom$indata.gene.mean <- NULL
  opossom$indata.original <- NULL
  opossom$indata.sample.mean <- NULL
  opossom$loglog.group.metadata <- NULL
  opossom$loglog.metadata <- NULL
  opossom$metadata <- NULL
  opossom$metagene.filter.list <- NULL
  opossom$n.0.m <- NULL
  opossom$output.paths <- NULL
  opossom$p.g.m <- NULL
  opossom$p.m <- NULL
  opossom$perc.DE.m <- NULL
  opossom$sd.g.m <- NULL
  opossom$som.nodes <- NULL
  opossom$som.result <- NULL
  opossom$secLvlSom.20.20 <- NULL
  opossom$secLvlSom.custom <- NULL
  opossom$t.g.m <- NULL
  opossom$t.m <- NULL
  opossom$groupwise.group.colors <- NULL
  opossom$unique.protein.ids <- NULL
  opossom$WAD.g.m <- NULL
  opossom$WAD.group.metadata <- NULL
  opossom$WAD.metadata <- NULL

  # Generate some additional letters
  opossom$LETTERS <- c(LETTERS, as.vector(sapply(1:10, function(x) {
    return(paste(LETTERS, x, sep=""))
  })))

  opossom$letters <- c(letters, as.vector(sapply(1:10, function(x) {
    return(paste(letters, x, sep=""))
  })))

  # Set default preferences
  opossom$preferences <- list(dataset.name = "Unnamed",
                              error.model = "all.samples",
                              dim.1stLvlSom = 20,
                              dim.2ndLvlSom = 20,
                              training.extension = 1,
                              rotate.SOM.portraits = 0,
                              flip.SOM.portraits = F,
                              database.dataset = "",
                              database.id.type = "",
                              geneset.analysis = F,
                              geneset.analysis.exact = F,
                              max.parallel.cores = detectCores() / 2,
                              spot.threshold.samples = 0.65,
                              spot.coresize.modules = 3,
                              spot.threshold.modules = 0.95,
                              spot.coresize.groupmap = 5,
                              spot.threshold.groupmap = 0.75,
                              feature.centralization = T,
                              sample.quantile.normalization = T,
                              pairwise.comparison.list = list())

  # Merge user supplied preferences
  opossom$preferences <-
    modifyList(opossom$preferences, preferences[names(opossom$preferences)])

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
  util.info("1SOM Dim:", opossom$preferences$dim.1stLvlSom)
  util.info("2SOM Dim:", opossom$preferences$dim.2ndLvlSom)

  # Prepare the environment
  if (!util.call(pipeline.prepare, opossom)) {
    return()
  }

  imagename <- paste(opossom$files.name, "pre.RData")
  util.info("Saving environment image:", imagename)
  save(opossom, file=imagename)

  # Execute the pipeline
  util.info("Processing Differential Expression")
  util.call(pipeline.calcStatistics, opossom)

  util.info("Detecting Spots")
  util.call(pipeline.detectSpotsSamples, opossom)
  util.call(pipeline.detectSpotsIntegral, opossom)
  util.call(pipeline.groupAssignment, opossom)

  util.info("Plotting Sample Portraits")
  util.call(pipeline.sampleExpressionPortraits, opossom)
  util.call(pipeline.sampleRankMaps, opossom)

  util.info("Processing Supporting Information")
  util.call(pipeline.supportingMaps, opossom)
  util.call(pipeline.entropyProfiles, opossom)
  util.call(pipeline.topologyProfiles, opossom)

  util.info("Processing 2nd level Metagene Analysis")
  dir.create(file.path(paste(opossom$files.name, "- Results"),
                       "2nd lvl Metagene Analysis"), showWarnings=F)

  util.call(pipeline.2ndLvlSimilarityAnalysis, opossom)
  util.call(pipeline.2ndLvlCorrelationAnalysis, opossom)
  util.call(pipeline.2ndLvlComponentAnalysis, opossom)
  util.call(pipeline.2ndLvlSom, opossom)


  if (opossom$preferences$geneset.analysis)
  {
    util.info("Processing Geneset Analysis")
    dir.create(paste(opossom$files.name, "- Results/Geneset Analysis"),
               showWarnings=F)

    util.call(pipeline.genesetStatisticSamples, opossom)
    util.call(pipeline.genesetStatisticIntegral, opossom)
    util.call(pipeline.genesetOverviews, opossom)

    util.info("Processing Geneset Profiles and Maps")
    util.call(pipeline.genesetProfilesAndMaps, opossom)

    util.info("Processing Cancer Hallmarks")
    util.call(pipeline.cancerHallmarks, opossom)

    util.info("Processing Chromosome Expression Reports")
    util.call(pipeline.chromosomeExpressionReports, opossom)
  }

  util.info("Processing Gene Lists")
  util.call(pipeline.geneLists, opossom)

  util.info("Processing Summary Sheets (Samples)")
  util.call(pipeline.summarySheetsSamples, opossom)

  util.info("Processing Summary Sheets (Spots)")
  util.call(pipeline.summarySheetsIntegral, opossom)

  util.info("Processing 3rd level Spot Analysis")
  dir.create(paste(opossom$files.name, "- Results/3rd lvl Spot Analysis"),
             showWarnings=F)

  util.call(pipeline.3rdLvlChromosomalEnrichment, opossom)
  util.call(pipeline.3rdLvlSummarySheets, opossom)
  util.call(pipeline.3rdLvlNetworks, opossom)

  util.info("Generating HTML Report")
  util.call(pipeline.htmlSummary, opossom)
  util.call(pipeline.htmlSampleSummary, opossom)
  util.call(pipeline.htmlIntegralSummary, opossom)
  util.call(pipeline.htmlGenesetAnalysis, opossom)

  # Save the opossom environment
  filename <- paste(opossom$files.name, ".RData", sep="")
  util.info("Saving environment image:", filename)
  save(opossom, file=filename)

  if (file.exists(imagename) && file.exists(filename))
  {
    file.remove(imagename)
  }

  # Run additional functions. (NOTE: They alter the environment)
  util.call(pipeline.groupAnalysis, opossom)
  util.call(pipeline.htmlGroupSummary, opossom)

  load(filename) # Reload opossom
  util.call(pipeline.differenceAnalyses, opossom)
  util.call(pipeline.htmlDifferencesSummary, opossom)

  load(filename) # Reload opossom
  util.call(pipeline.signatureSets, opossom)

  util.info("Finished:", format(Sys.time(), "%a %b %d %X"))
}
