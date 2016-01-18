# Creates a new opossom environment
opossom.new <- function(preferences=NULL)
{
  # Init the environment
  env <- new.env()
  env$t.ensID.m <- NULL
  env$colramp <- NULL
  env$Fdr.g.m <- NULL
  env$fdr.g.m <- NULL
  env$files.name <- NULL
  env$gene.descriptions <- NULL
  env$gene.ids <- NULL
  env$gene.names <- NULL
  env$gene.positions <- NULL
  env$gene.positions.list <- NULL
  env$gene.positions.table <- NULL
  env$gene.coordinates <- NULL
  env$group.silhouette.coef <- NULL
  env$group.colors <- NULL
  env$group.labels <- NULL
  env$gs.def.list <- NULL
  env$samples.GSZ.scores <- NULL
  env$spot.list.correlation <- NULL
  env$spot.list.dmap <- NULL
  env$spot.list.group.overexpression <- NULL
  env$spot.list.kmeans <- NULL
  env$spot.list.overexpression <- NULL
  env$spot.list.samples <- NULL
  env$spot.list.underexpression <- NULL
  env$indata <- NULL
  env$indata.gene.mean <- NULL
  env$indata.sample.mean <- NULL
  env$metadata <- NULL
  env$metagene.filter.list <- NULL
  env$n.0.m <- NULL
  env$output.paths <- NULL
  env$pat.labels <- NULL
  env$p.g.m <- NULL
  env$p.m <- NULL
  env$perc.DE.m <- NULL
  env$som.nodes <- NULL
  env$som.result <- NULL
  env$secLvlSom.20.20 <- NULL
  env$secLvlSom.custom <- NULL
  env$t.g.m <- NULL
  env$t.m <- NULL
  env$groupwise.group.colors <- NULL
  env$unique.protein.ids <- NULL
  env$WAD.g.m <- NULL

  # Generate some additional letters
  env$LETTERS <- c(LETTERS, as.vector(sapply(1:10, function(x) {
    return(paste(LETTERS, x, sep=""))
  })))

  env$letters <- c(letters, as.vector(sapply(1:10, function(x) {
    return(paste(letters, x, sep=""))
  })))

  # Set default preferences
  env$preferences <- list(dataset.name = "Unnamed",
                          dim.1stLvlSom = "auto",
                          dim.2ndLvlSom = 20,
                          training.extension = 1,
                          rotate.SOM.portraits = 0,
                          flip.SOM.portraits = FALSE,
                          database.dataset = "auto",
                          database.id.type = "",
                          geneset.analysis = TRUE,
                          geneset.analysis.exact = FALSE,
                          spot.coresize.modules = 3,
                          spot.threshold.modules = 0.95,
                          spot.coresize.groupmap = 5,
                          spot.threshold.groupmap = 0.75,
                          feature.centralization = TRUE,
                          sample.quantile.normalization = TRUE,
                          pairwise.comparison.list = list())

  # Merge user supplied information
  if (!is.null(preferences))
  {
    env$preferences <-
      modifyList(env$preferences, preferences[names(env$preferences)])
  }
  if(!is.null(preferences$indata))
  {
    env$indata <- preferences$indata
  }
  if(!is.null(preferences$group.labels))
  {
    env$group.labels <- preferences$group.labels
  }
  if(!is.null(preferences$group.colors))
  {
    env$group.colors <- preferences$group.colors
  }
  
  return(env)
}

# Executes the oposSOM pipeline.
opossom.run <- function(env)
{
  env$preferences$system.info <- Sys.info()
  env$preferences$session.info <- sessionInfo()
  env$preferences$started <- format(Sys.time(), "%a %d %b %Y %X")

  # Output some info
  util.info("Started:", env$preferences$started)
  util.info("Name:", env$preferences$dataset.name)

  # Dump frames on error
  error.option <- getOption("error")
  options(error=quote({dump.frames(to.file=TRUE)}))

  # Prepare the environment
  if (!util.call(pipeline.prepare, env)) {
    return()
  }

  imagename <- paste(env$files.name, "pre.RData")
  util.info("Saving environment image:", imagename)
  save(env, file=imagename)

  # Execute the pipeline
  util.info("Processing Differential Expression")
  util.call(pipeline.calcStatistics, env)

  util.info("Detecting Spots")
  util.call(pipeline.detectSpotsSamples, env)
  util.call(pipeline.detectSpotsIntegral, env)
  util.call(pipeline.patAssignment, env)
  util.call(pipeline.groupAssignment, env)

  util.info("Plotting Sample Portraits")
  util.call(pipeline.sampleExpressionPortraits, env)
  util.call(pipeline.sampleRankMaps, env)

  util.info("Processing Supporting Information")
  util.call(pipeline.supportingMaps, env)
  util.call(pipeline.entropyProfiles, env)
  util.call(pipeline.topologyProfiles, env)

  util.info("Processing 2nd level Metagene Analysis")
  dir.create(file.path(paste(env$files.name, "- Results"),
                       "2nd lvl Metagene Analysis"), showWarnings=FALSE)

  util.call(pipeline.2ndLvlSimilarityAnalysis, env)
  util.call(pipeline.2ndLvlCorrelationAnalysis, env)
  util.call(pipeline.2ndLvlComponentAnalysis, env)
  util.call(pipeline.2ndLvlSom, env)


  if (env$preferences$geneset.analysis)
  {
    util.info("Processing Geneset Analysis")
    dir.create(paste(env$files.name, "- Results/Geneset Analysis"),
               showWarnings=FALSE)

    util.call(pipeline.genesetStatisticSamples, env)
    util.call(pipeline.genesetStatisticIntegral, env)
    util.call(pipeline.genesetOverviews, env)

    util.info("Processing Geneset Profiles and Maps")
    util.call(pipeline.genesetProfilesAndMaps, env)

    util.info("Processing Cancer Hallmarks")
    util.call(pipeline.cancerHallmarks, env)

    util.info("Processing Chromosome Expression Reports")
    util.call(pipeline.chromosomeExpressionReports, env)
  }

  util.info("Processing Gene Lists")
  util.call(pipeline.geneLists, env)

  util.info("Processing Summary Sheets (Samples)")
  util.call(pipeline.summarySheetsSamples, env)

  util.info("Processing Summary Sheets (Spots)")
  util.call(pipeline.summarySheetsIntegral, env)

  util.info("Processing 3rd level Spot Analysis")
  dir.create(paste(env$files.name, "- Results/3rd lvl Spot Analysis"),
             showWarnings=FALSE)

  util.call(pipeline.3rdLvlChromosomalEnrichment, env)
  util.call(pipeline.3rdLvlSummarySheets, env)
  util.call(pipeline.3rdLvlNetworks, env)

  util.info("Generating HTML Report")
  util.call(pipeline.htmlSummary, env)
  util.call(pipeline.htmlSampleSummary, env)
  util.call(pipeline.htmlIntegralSummary, env)
  util.call(pipeline.htmlGenesetAnalysis, env)

  # Save the opossom environment
  filename <- paste(env$files.name, ".RData", sep="")
  util.info("Saving environment image:", filename)
  save(env, file=filename)

  if (file.exists(imagename) && file.exists(filename))
  {
    file.remove(imagename)
  }

  # Run additional functions. (NOTE: They alter the environment)
  util.call(pipeline.summarySheetsPATs, env)
  
  load(filename) # Reload env
  util.call(pipeline.groupAnalysis, env)
  util.call(pipeline.htmlGroupSummary, env)

  load(filename) # Reload env
  util.call(pipeline.differenceAnalyses, env)
  util.call(pipeline.htmlDifferencesSummary, env)

  # Restore old error behaviour
  options(error=error.option)

  util.info("Finished:", format(Sys.time(), "%a %b %d %X"))
}
