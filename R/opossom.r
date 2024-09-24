# Creates a new opossom environment
opossom.new <- function(preferences=NULL)
{
  # Init the environment
  env <- new.env()
  env$color.palette.portraits <- NULL
  env$color.palette.heatmaps <- NULL
  env$fdr.g.m <- NULL
  env$files.name <- NULL
  env$gene.info <- NULL
  env$chromosome.list <- NULL
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
  env$spot.list.underexpression <- NULL
  env$indata <- NULL
	env$indata.ensID.m <- NULL
  env$indata.gene.mean <- NULL
  env$indata.sample.mean <- NULL
  env$metadata <- NULL
  env$n.0.m <- NULL
  env$output.paths <- NULL
  env$pat.labels <- NULL
	env$pat.colors <- NULL
  env$p.g.m <- NULL
  env$p.m <- NULL
  env$perc.DE.m <- NULL
  env$psf.results.samples <- NULL
  env$psf.results.groups <- NULL
  env$som.result <- NULL
  env$groupwise.group.colors <- NULL
  env$csv.function <- write.csv2

  # Generate some additional letters
  env$LETTERS <- c(LETTERS, as.vector(sapply(1:10, function(x) {
    return(paste(LETTERS, x, sep=""))
  })))

  env$letters <- c(letters, as.vector(sapply(1:10, function(x) {
    return(paste(letters, x, sep=""))
  })))

  # Set default preferences
  env$preferences <- list(dataset.name = "Unnamed",
													note = "",
													max.cores = detectCores()-1,
                          dim.1stLvlSom = "auto",
                          dim.2ndLvlSom = 20,
                          training.extension = 1,
                          rotate.SOM.portraits = 0,
                          flip.SOM.portraits = FALSE,
													colorblindsave.portraits = TRUE,
                          activated.modules = list( "reporting" = TRUE,
                                                    "primary.analysis" = TRUE, 
                                                    "sample.similarity.analysis" = TRUE,
                                                    "geneset.analysis" = TRUE, 
                                                    "psf.analysis" = TRUE,
                                                    "group.analysis" = TRUE,
                                                    "difference.analysis" = TRUE,
                                                    "largedata.mode" = NULL),
                          database.biomart = "ENSEMBL_MART_ENSEMBL",
                          database.host = "https://jan2020.archive.ensembl.org",
                          database.dataset = "auto",
                          database.id.type = "",
                          standard.spot.modules = "dmap",
                          spot.coresize.modules = 3,
                          spot.threshold.modules = 0.95,
                          spot.coresize.groupmap = 5,
                          spot.threshold.groupmap = 0.75,
                          adjust.autogroup.number = 0,
                          feature.centralization = TRUE,
                          sample.quantile.normalization = TRUE,
                          pairwise.comparison.list = NULL)

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
  util.info("oposSOM is ready to fly! Starting analysis.")
  util.info("Name:", env$preferences$dataset.name)

  #### Preparation & Calculation part ####
  env <- pipeline.checkInputParameters(env)
  if (!env$passedInputChecking) {
    return()
  }
  
  if(env$preferences$activated.modules$primary.analysis)
  {
    env$preferences$system.info <- Sys.info()
    env$preferences$session.info <- sessionInfo()
    env$preferences$started <- format(Sys.time(), "%a %d %b %Y %X")
  }
  
  if(env$preferences$activated.modules$reporting)
  {
    # create output dirs
    dir.create(paste(env$files.name, "- Results"), showWarnings=FALSE)
    dir.create(paste(env$files.name, "- Results/CSV Sheets"), showWarnings=FALSE)
  	setwd(paste(env$files.name, "- Results"))
	
    if(env$preferences$activated.modules$primary.analysis)
    {
      pipeline.qualityCheck(env)
    } 
  }
  if(env$preferences$activated.modules$primary.analysis || env$preferences$activated.modules$geneset.analysis)
  {
    util.info("Loading gene annotation data.")
    env <- pipeline.prepareAnnotation(env)
  }
  
  if(env$preferences$activated.modules$primary.analysis)
  {
    util.info("Processing SOM. This may take several time until next notification.")
    env <- pipeline.prepareIndata(env)
    env <- pipeline.generateSOM(env)
    
    filename <- paste(env$files.name, "pre.RData")
    util.info("Saving environment image:", filename)
    save(env, file=filename)
    
    if( !env$preferences$activated.modules$largedata.mode )
    {
      util.info("Processing Differential Expression Statistics")
      env <- pipeline.diffExpressionStatistics(env)
    }

    util.info("Detecting Spots")
    env <- pipeline.detectSpotsModules(env)
    env <- pipeline.patAssignment(env)
    env <- pipeline.groupAssignment(env)
  }

  if (env$preferences$activated.modules$geneset.analysis)
  {
    util.info("Calculating Geneset Enrichment")
    env <- pipeline.genesetStatisticSamples(env)
    env <- pipeline.genesetStatisticModules(env)
  }
  
  if (env$preferences$activated.modules$psf.analysis)
  {
    util.info("Calculating Pathway Signal Flow (PSF)")
    env <- pipeline.PSFcalculation(env)    
  }
  
  if(env$preferences$activated.modules$primary.analysis || env$preferences$activated.modules$geneset.analysis)
  {    
    filename <- paste(env$files.name, ".RData", sep="")
    util.info("Saving environment image:", filename)
    save(env, file=filename)
    
    if (file.exists(paste(env$files.name, "pre.RData")) && file.exists(filename))
    {
      file.remove(paste(env$files.name, "pre.RData"))
    }
  }  
    
  #### Reporting part ####
  
  if(env$preferences$activated.modules$reporting)
  {
  
    util.info("Plotting Supporting Information")
    pipeline.supportingMaps(env)
    pipeline.entropyProfiles(env)
    pipeline.topologyProfiles(env)

    
    if(length(env$chromosome.list) > 0)
    {
      util.info("Plotting Chromosome Expression Reports")
      pipeline.chromosomeExpressionReports(env)
    }
    
    if( !env$preferences$activated.modules$largedata.mode && ncol(env$indata)<1000 )
    {
      util.info("Plotting Sample Portraits")
      pipeline.sampleExpressionPortraits(env)
    } 
    
    if ( env$preferences$activated.modules$sample.similarity.analysis && 
         ncol(env$indata) > 2 &&
         !env$preferences$activated.modules$largedata.mode )
    {    
      util.info("Plotting Sample Similarity Analysis")
      dir.create("Sample Similarity Analysis", showWarnings=FALSE)
      
      pipeline.sampleSimilarityAnalysisED(env)
      pipeline.sampleSimilarityAnalysisCor(env)
      pipeline.sampleSimilarityAnalysisICA(env)
      pipeline.sampleSimilarityAnalysisSOM(env)
    }
    
    if (env$preferences$activated.modules$geneset.analysis &&
        !env$preferences$activated.modules$largedata.mode )
    {
      dir.create("Geneset Analysis", showWarnings=FALSE)
      
      util.info("Plotting Geneset Enrichment Heatmaps")
      pipeline.genesetOverviews(env)
      
      util.info("Plotting Geneset Profiles and Maps")
      pipeline.genesetProfilesAndMaps(env)
      
      util.info("Calculating Cancer Hallmark Enrichment")
      pipeline.cancerHallmarks(env)
    }
    
    if (env$preferences$activated.modules$psf.analysis)
    {
      util.info("Plotting PSF results")
      pipeline.PSFoutput(env)
    }
    
    util.info("Writing Gene Lists")
    pipeline.geneLists(env)

    if( !env$preferences$activated.modules$largedata.mode && ncol(env$indata)<1000 )
    {
      util.info("Plotting Summary Sheets (Samples)")
      pipeline.summarySheetsSamples(env)
    }
    
    util.info("Plotting Summary Sheets (Modules & PATs)")
    pipeline.summarySheetsModules(env)
    pipeline.summarySheetsPATs(env)
      

    if(env$preferences$activated.modules$group.analysis && length(unique(env$group.labels)) >= 2)
    {
      util.info("Processing Group-centered Analyses")
      pipeline.groupAnalysis(env)
    }
  
    if(env$preferences$activated.modules$difference.analysis)
    {
      util.info("Processing Difference Analyses")
      pipeline.differenceAnalyses(env)
    }

    util.info("Generating HTML Report")
    pipeline.htmlSampleSummary(env)
    pipeline.htmlModuleSummary(env)
    pipeline.htmlGenesetAnalysis(env)  
    pipeline.htmlPsfAnalysis(env)
    pipeline.htmlSummary(env)
    
  }    
    
  util.info("Finished:", format(Sys.time(), "%a %b %d %X"))
	
	return(env)
}
