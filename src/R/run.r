pipeline.run <- function() {

  # output some info
  util.info("Started:", format(Sys.time(), "%a %b %d %X"))
  util.info("Setting:", preferences$dataset.name)
  util.info("1SOM Dim:", preferences$dim.som1)
  util.info("2SOM Dim:", preferences$dim.som2)

  # check input parameters/data

  if (is.null(indata))
  {
    util.fatal("No indata supplied!")
    return()
  }

  if (class(indata) != "matrix" || mode(indata) != "numeric" || storage.mode(indata) != "numeric")
  {
    rn = rownames(indata)
    indata <<- apply(indata, 2, function(x){ as.numeric(as.vector(x)) })
    rownames(indata) <<- rn
    storage.mode(indata) <<- "numeric"
    #util.warn("Converted indata to numerical matrix")
  }

  if (length(rownames(indata)) == 0)
  {
    rownames(indata) <<- as.character(1:nrow(indata))
    preferences$geneset.analysis <<- F
    util.warn("No rownames found. Set them to 1,2,3,4...")
  }

  if (length(colnames(indata)) == 0)
  {
    colnames(indata) <<- paste("Sample", c(1:ncol(indata)))
    util.warn("No colnames found. Set them to 1,2,3,4...")
  }

  if (any(duplicated(rownames(indata))))
  {
    indata <<- do.call(rbind, by(indata, rownames(indata), colMeans))[unique(rownames(indata)),]
    util.warn("Duplicate rownames. Averaged multiple features")
  }

  na.rows <- which(apply(apply(indata, 1, is.na), 2, sum) > 0)

  if (length(na.rows) > 0)
  {
    indata <<- indata[-na.rows,]

    if (exists("indata.original")) {
      indata.original <<- indata.original[-na.rows,]
    }
    util.warn("Removed NAs from data set")
  }

  ## set up global variables

  files.name <<- preferences$dataset.name

  while (file.exists(paste(files.name, ".RData", sep=""))) {
    files.name <<- paste(files.name, "+", sep="")
  }

  output.paths <<- c("LPE" = paste(files.name, "- Results/LPE"),
                     "CSV" = paste(files.name, "- Results/CSV Sheets"),
                     "Summary Sheets Samples"= paste(files.name, "- Results/Summary Sheets - Samples"),
                     "Summary Sheets Integral"= paste(files.name, "- Results/Summary Sheets - Integral"))

  # create output dirs
  dir.create(paste(files.name, "- Results"), showWarnings=F)
  dir.create(paste(files.name, "- Results/CSV Sheets"), showWarnings=F)

  if (is.null(colramp)) {
    colramp <<- colorRampPalette(c("darkblue", "blue", "lightblue", "green",
                                   "yellow", "red", "darkred"))
  }

  if (length(group.labels) != ncol(indata) || length(group.colors) != ncol(indata))
  {
    group.labels <<- NULL
    group.colors <<- NULL
    util.warn("Group assignment doesnt fit number of samples")
  }

  if (!is.null(group.labels) && max(table(group.labels)) == 1)
  {
    group.labels <<- NULL
    group.colors <<- NULL
    util.warn("Each sample has an own group")
  }

  if (!is.null(group.labels))
  {
    group.labels <<- as.character(group.labels)
    names(group.labels) <<- colnames(indata)

    if (is.null(group.colors))
    {
      group.colors <<- rep("", ncol(indata))

      for (i in 1:length(unique(group.labels)))
      {
        group.colors[which(group.labels == unique(group.labels)[i])] <<-
          colorRampPalette(c("blue3", "blue", "lightblue", "green", "gold", "red", "red3"))(length(unique(group.labels)))[i]
      }
    }

    # catch userdefined group.colors --> convert to #hex
    if (length(unique(substr(group.colors, 1, 1)) > 1) || unique(substr(group.colors, 1, 1))[1] != "#")
    {
      group.colors <<- apply(col2rgb(group.colors), 2, function(x) { rgb(x[1]/255, x[2]/255, x[3]/255) })
    }
    names(group.colors) <<- colnames( indata )
  } else
  {
    group.labels <<- rep("sample", ncol(indata))
    names(group.labels) <<- colnames(indata)

    group.colors <<- colramp(ncol(indata))
    names(group.colors) <<- colnames(indata)
  }

  unique.group.colors <<- group.colors[match(unique(group.labels), group.labels)]
  names(unique.group.colors) <<- unique(group.labels)

  # prepare data

  indata.sample.mean <<- colMeans(indata)
  environment(pipeline.qualityCheck) <- environment()
  pipeline.qualityCheck()

  if (preferences$sample.quantile.normalization)
  {
    indata <<- Quantile.Normalization(indata)

    if (!is.null(indata.original))
    {
      indata.original <<- Quantile.Normalization(indata.original)
      util.warn("Separate quantile normalization of indata AND indata.original")
    }
  }

  if (!is.null(indata.original) && any(dim(indata) != dim(indata.original)))
  {
    indata.original <<- NULL
    util.warn("Existing 'indata.original' does not fit 'indata' object")
  }

  if (is.null(indata.original))
  {
    indata.original <<- indata
  }

  if (preferences$error.model == "replicates")
  {
    indata <<- do.call(cbind, by(t(indata), colnames(indata), mean))[,unique(colnames(indata))]
    group.labels <<- group.labels[colnames(indata)]
    group.colors <<- group.colors[colnames(indata)]
    indata.sample.mean <<- tapply(indata.sample.mean, colnames(indata.original), mean)[colnames(indata)]
  } else
  {
    colnames(indata) <<- make.unique(colnames(indata))
    names(group.labels) <<- make.unique(names(group.labels))
    names(group.colors) <<- make.unique(names(group.colors))
  }

  indata.mean.level <<- rowMeans(indata)

  if (preferences$feature.mean.normalization)
  {
    indata <<- indata - indata.mean.level
  }

  util.info("Load Annotation Data")
  environment(pipeline.prepareAnnotation) <- environment()
  pipeline.prepareAnnotation()


  ## SOM
  util.info("Processing SOM")

  som.result <<- som.init(indata, xdim=preferences$dim.som1, ydim=preferences$dim.som1, init="linear")

  # Rotate/Flip First lvl SOMs

  if (preferences$rotate.som1 > 0)
  {
    for(i in 1:preferences$rotate.som1)
    {
      o <- matrix(c(1:(preferences$dim.som1^2)), preferences$dim.som1, preferences$dim.som1, byrow=T)
      o <- o[rev(1:preferences$dim.som1),]
      som.result <<- som.result[as.vector(o),]
    }
  }

  if(preferences$flip.som1)
  {
    o <- matrix(c(1:(preferences$dim.som1^2)), preferences$dim.som1, preferences$dim.som1, byrow=T)
    som.result <<- som.result[as.vector(o),]
  }


  # Train SOM

  t1 <- system.time({
    som.result <<- som.train(indata, som.result, xdim=preferences$dim.som1,
                             ydim=preferences$dim.som1, alpha=0.05,
                             radius=preferences$dim.som1,
                             rlen=nrow(indata)*2*preferences$training.extension,
                             inv.alp.c=nrow(indata)*2*preferences$training.extension/100)
  })
  util.info("Remaining ~", ceiling(5*t1[3]/60), "min ~", round(5*t1[3]/3600,1),"h")

  som.result <<- som.train(indata, som.result$code, xdim=preferences$dim.som1,
                           ydim=preferences$dim.som1, alpha=0.02,
                           radius=min(3, preferences$dim.som1),
                           rlen=nrow(indata)*10*preferences$training.extension,
                           inv.alp.c=nrow(indata)*10*preferences$training.extension/100)

  # TODO Can we throw this in the bin?
  #som.result <<- som( indata, xdim=preferences$dim.som1, ydim=preferences$dim.som1 )

  metadata <<- som.result$code
  colnames(metadata) <<- colnames(indata)

  som.result$code <<- NA

  loglog.metadata <<- apply(metadata, 2, function(x)
  {
    meta.sign <- sign(x)
    meta <- log10(abs(x))
    meta <- meta - min( meta, na.rm=T )
    return(meta * meta.sign)
  })

  WAD.metadata <<- apply(metadata ,2, function(x) { x * ((x - min(x)) / (max(x) - min(x))) })


  ##  Group SOMs

  if (length(unique(group.labels)) > 1) # mean group metagenes
  {
    group.metadata <<- do.call(cbind, by(t(metadata), group.labels, colMeans))[,unique(group.labels)]

    loglog.group.metadata <<- do.call(cbind, by(t(loglog.metadata), group.labels, colMeans))[,unique(group.labels)]

    WAD.group.metadata <<- do.call(cbind, by(t(WAD.metadata), group.labels, colMeans))[,unique(group.labels)]
  }


  ## set up SOM dependent variables
  genes.coordinates <<- apply(som.result$visual[,c(1,2)]+1, 1, paste, collapse=" x ")
  names(genes.coordinates) <<- rownames(indata)

  som.nodes <<- (som.result$visual[,"x"] + 1) + som.result$visual[,"y"] * preferences$dim.som1
  names(som.nodes) <<- rownames(indata)


  ## Output

  util.info("Saving workspace image:", files.name, "pre.RData")
  # TODO Save opossom environment only?
  save.image(paste(files.name, " pre.RData" , sep=""))

  util.info("Processing Differential Expression")
  environment(pipeline.calcStatistics) <- environment()
  pipeline.calcStatistics()

  util.info("Detecting Spots")
  environment(pipeline.detectSpotsSamples) <- environment()
  pipeline.detectSpotsSamples()

  environment(pipeline.detectSpotsIntegral) <- environment()
  pipeline.detectSpotsIntegral()

############### TODO ###

  source("R/source/group_assignment.r", local=TRUE)


  cat( "Plotting Sample Portraits\n" ); flush.console()
  source("R/source/sample_expression_portraits.r", local=TRUE)
  source("R/source/sample_rank_maps.r", local=TRUE)


  cat( "Processing Supporting Information\n" ); flush.console()
  source("R/source/supporting_maps.r", local=TRUE)
  source("R/source/entropy_profiles.r", local=TRUE)
  source("R/source/topology_profiles.r", local=TRUE)


  cat( "Processing 2nd level Metagene Analysis\n" ); flush.console()
  dir.create( paste( files.name, "- Results/2nd lvl Metagene Analysis" ), showWarnings=F )

  source("R/source/2nd_lvl_similarity_analysis.r", local=TRUE)
  source("R/source/2nd_lvl_correlation_analysis.r", local=TRUE)
  source("R/source/2nd_lvl_component_analysis.r", local=TRUE)
  source("R/source/2nd_lvl_som.r", local=TRUE)



  if( preferences$geneset.analysis )
  {
    dir.create( paste( files.name, "- Results/Geneset Analysis" ), showWarnings=F )

    source("R/source/geneset_statistic_samples.r", local=TRUE)
    source("R/source/geneset_statistic_integral.r", local=TRUE)
    source("R/source/geneset_overviews.r", local=TRUE)
    source("R/source/geneset_profiles_and_maps.r", local=TRUE)

    source("R/source/cancer_hallmarks.r", local=TRUE)
    source("R/source/chromosome_expression_reports.r", local=TRUE)
  }



  cat( "Gene Lists\n" ); flush.console()
  source("R/source/gene_lists.r", local=TRUE)

  cat( "Summary Sheets: Samples\n" ); flush.console()
  source("R/source/summary_sheets_samples.r", local=TRUE)

  cat( "Summary Sheets: Spots\n" ); flush.console()
  source("R/source/summary_sheets_integral.r", local=TRUE)



  cat( "Processing 3rd level Spot Analysis\n" ); flush.console()
  dir.create( paste( files.name, "- Results/3rd lvl Spot Analysis" ), showWarnings=F )

  source("R/source/3rd_lvl_chromosomal_enrichment.r", local=TRUE)
  source("R/source/3rd_lvl_summary_sheets.r", local=TRUE)
  source("R/source/3rd_lvl_networks.r", local=TRUE)




  cat( "Generating HTML Report\n" ); flush.console()
  source("R/source/html_summary.r", local=TRUE)
  source("R/source/html_sample_summary.r", local=TRUE)
  source("R/source/html_integral_summary.r", local=TRUE)
  source("R/source/html_geneset_analysis.r", local=TRUE)



  cat( "Clean and store Workspace\n" ); flush.console()
  source("R/source/workspace_cleanup.r", local=TRUE)
  save.image( paste( files.name, ".RData" , sep="" ) )


  if( file.exists( paste( files.name, " pre.RData" , sep="") ) && file.exists( paste( files.name, ".RData" , sep="") ) )
    r = file.remove( paste( files.name, " pre.RData" , sep="") )




  ## additional scripts

  source("R/source/group_analyses.r", local=TRUE)
  source("R/source/difference_analyses.r", local=TRUE)



#  cat( "Spot Filtering\n" ); flush.console()
#  source("R/source/3rd_lvl_overexpression_genenet.r", local=TRUE)
#  source("R/source/3rd_lvl_spot_filter.r", local=TRUE)
  source("R/source/signature_sets.r", local=TRUE)



  cat( "Finished:", format(Sys.time(), "%a %b %d %X\n\n" ) )
  flush.console()





}
