pipeline.prepare <- function()
{
  ## check preferences
  if (!is.character(preferences$dataset.name))
  {
    util.warn("Invalid value of \"dataset.name\". Using \"Unnamed\"")
    preferences$dataset.name <<- "Unnamed"
  }

  if (!preferences$error.model %in% c("replicates", "all.samples", "groups")) {
    util.warn("Invalid value of \"error.model\". Using \"all.samples\"")
    preferences$error.model <<- "all.samples"
  }

  if (!is.numeric(preferences$dim.1stLvlSom) || preferences$dim.1stLvlSom < 1)
  {
    util.warn("Invalid value of \"dim.1stLvlSom\". Using 20")
    preferences$dim.1stLvlSom <<- 20
  }

  if (!is.numeric(preferences$dim.2ndLvlSom) || preferences$dim.2ndLvlSom < 1)
  {
    util.warn("Invalid value of \"dim.2ndLvlSom\". Using 20")
    preferences$dim.2ndLvlSom <<- 20
  }

  if (!is.numeric(preferences$training.extension) ||
      preferences$training.extension < 1 ||
      preferences$training.extension > 10)
  {
    util.warn("Invalid value of \"training.extension\". Using 1")
    preferences$training.extension <<- 1
  }

  if (!is.numeric(preferences$rotate.som1) ||
      preferences$rotate.som1 < 0 ||
      preferences$rotate.som1 > 4)
  {
    util.warn("Invalid value of \"rotate.som1\". Using 0")
    preferences$rotate.som1 <<- 0
  }

  if (!is.logical(preferences$flip.SOM.portraits))
  {
    util.warn("Invalid value of \"flip.SOM.portraits\". Using FALSE")
    preferences$flip.SOM.portraits <<- F
  }

  if (!is.character(preferences$database.dataset))
  {
    util.warn("Invalid value of \"database.dataset\". Using \"\"")
    preferences$database.dataset <<- ""
  }

  if (!is.character(preferences$database.id.type))
  {
    util.warn("Invalid value of \"database.id.type\". Using \"\"")
    preferences$database.id.type <<- ""
  }

  if (!is.logical(preferences$geneset.analysis))
  {
    util.warn("Invalid value of \"geneset.analysis\". Using FALSE")
    preferences$geneset.analysis <<- F
  }

  if (!is.logical(preferences$geneset.analysis.exact))
  {
    util.warn("Invalid value of \"geneset.analysis.exact\". Using FALSE")
    preferences$geneset.analysis.exact <<- F
  }

  if (!is.numeric(preferences$max.parallel.cores) ||
      preferences$max.parallel.cores < 1)
  {
    preferences$max.parallel.cores <<- detectCores() / 2
    util.warn("Invalid value of \"max.parallel.cores\". Using",
              preferences$max.parallel.cores)
  }

  if (!is.numeric(preferences$sample.spot.cutoff) ||
      preferences$sample.spot.cutoff <= 0 ||
      preferences$sample.spot.cutoff >= 1)
  {
    util.warn("Invalid value of \"sample.spot.cutoff\". Using 0.65")
    preferences$sample.spot.cutoff <<- 0.65
  }

  if (!is.numeric(preferences$summary.spot.core) ||
      preferences$summary.spot.core < 1 ||
      preferences$summary.spot.core > 10)
  {
    util.warn("Invalid value of \"summary.spot.core\". Using 3")
    preferences$summary.spot.core <<- 3
  }

  if (!is.numeric(preferences$summary.spot.threshold) ||
      preferences$summary.spot.threshold <= 0 ||
      preferences$summary.spot.threshold >= 1)
  {
    util.warn("Invalid value of \"summary.spot.threshold\". Using 0.95")
    preferences$summary.spot.threshold <<- 0.95
  }

  if (!is.numeric(preferences$group.spot.core) ||
      preferences$group.spot.core < 1 ||
      preferences$group.spot.core > 10)
  {
    util.warn("Invalid value of \"group.spot.core\". Using 5")
    preferences$group.spot.core <<- 5
  }

  if (!is.numeric(preferences$group.spot.threshold) ||
      preferences$group.spot.threshold <= 0 ||
      preferences$group.spot.threshold >= 1)
  {
    util.warn("Invalid value of \"group.spot.threshold\". Using 0.75")
    preferences$group.spot.threshold <<- 0.75
  }

  if (!is.logical(preferences$feature.centralization))
  {
    util.warn("Invalid value of \"feature.centralization\". Using TRUE")
    preferences$feature.centralization <<- T
  }

  if (!is.logical(preferences$sample.quantile.normalization))
  {
    util.warn("Invalid value of \"sample.quantile.normalization\". Using TRUE")
    preferences$sample.quantile.normalization <<- T
  }

  # check input parameters/data
  if (is.null(indata))
  {
    util.fatal("No indata supplied!")
    return(F)
  }

  if (class(indata) != "matrix" && (is.null(dim(indata)) || dim(indata) < 1))
  {
    util.fatal("Invalid indata supplied!")
    return(F)
  }

  if (class(indata) != "matrix" || mode(indata) != "numeric" || storage.mode(indata) != "numeric")
  {
    rn <- rownames(indata)
    indata <<- apply(indata, 2, function(x){ as.numeric(as.vector(x)) })
    rownames(indata) <<- rn
    storage.mode(indata) <<- "numeric"
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

    if (!is.null(indata.original)) {
      indata.original <<- indata.original[-na.rows,]
    }
    util.warn("Removed NAs from data set")
  }

  ## set up global variables

  files.name <<- preferences$dataset.name

  while (file.exists(paste(files.name, ".RData", sep=""))) {
    files.name <<- paste(files.name, "+", sep="")
  }

  output.paths <<-
    c("LPE"=paste(files.name, "- Results/LPE"),
      "CSV"=paste(files.name, "- Results/CSV Sheets"),
      "Summary Sheets Samples"=paste(files.name, "- Results/Summary Sheets - Samples"),
      "Summary Sheets Integral"=paste(files.name, "- Results/Summary Sheets - Integral"))

  # create output dirs
  dir.create(paste(files.name, "- Results"), showWarnings=F)
  dir.create(paste(files.name, "- Results/CSV Sheets"), showWarnings=F)

  if (is.null(colramp)) {
    colramp <<- colorRampPalette(c("darkblue", "blue", "lightblue", "green",
                                   "yellow", "red", "darkred"))
  }

  if ((!is.null(group.labels) && length(group.labels) != ncol(indata)) ||
      (!is.null(group.colors) && length(group.colors) != ncol(indata)))
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
    names(group.colors) <<- colnames(indata)
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

  if (preferences$feature.centralization)
  {
    indata <<- indata - indata.mean.level
  }

  util.info("Load Annotation Data")
  environment(pipeline.prepareAnnotation) <- environment()
  pipeline.prepareAnnotation()


  ## SOM
  util.info("Processing SOM")

  som.result <<- som.init(indata, xdim=preferences$dim.1stLvlSom, ydim=preferences$dim.1stLvlSom, init="linear")

  # Rotate/Flip First lvl SOMs

  if (preferences$rotate.som1 > 0)
  {
    for (i in 1:preferences$rotate.som1)
    {
      o <- matrix(c(1:(preferences$dim.1stLvlSom^2)), preferences$dim.1stLvlSom, preferences$dim.1stLvlSom, byrow=T)
      o <- o[rev(1:preferences$dim.1stLvlSom),]
      som.result <<- som.result[as.vector(o),]
    }
  }

  if (preferences$flip.SOM.portraits)
  {
    o <- matrix(c(1:(preferences$dim.1stLvlSom^2)), preferences$dim.1stLvlSom, preferences$dim.1stLvlSom, byrow=T)
    som.result <<- som.result[as.vector(o),]
  }


  # Train SOM

  t1 <- system.time({
    som.result <<- som.train(indata, som.result, xdim=preferences$dim.1stLvlSom,
                             ydim=preferences$dim.1stLvlSom, alpha=0.05,
                             radius=preferences$dim.1stLvlSom,
                             rlen=nrow(indata)*2*preferences$training.extension,
                             inv.alp.c=nrow(indata)*2*preferences$training.extension/100)
  })
  util.info("Remaining ~", ceiling(5*t1[3]/60), "min ~", round(5*t1[3]/3600,1),"h")

  som.result <<- som.train(indata, som.result$code, xdim=preferences$dim.1stLvlSom,
                           ydim=preferences$dim.1stLvlSom, alpha=0.02,
                           radius=min(3, preferences$dim.1stLvlSom),
                           rlen=nrow(indata)*10*preferences$training.extension,
                           inv.alp.c=nrow(indata)*10*preferences$training.extension/100)

  # TODO Can we throw this in the bin?
  #som.result <<- som(indata, xdim=preferences$dim.1stLvlSom, ydim=preferences$dim.1stLvlSom)

  metadata <<- som.result$code
  colnames(metadata) <<- colnames(indata)

  som.result$code <<- NA

  loglog.metadata <<- apply(metadata, 2, function(x)
  {
    meta.sign <- sign(x)
    meta <- log10(abs(x))
    meta <- meta - min(meta, na.rm=T)
    return(meta * meta.sign)
  })

  WAD.metadata <<- apply(metadata ,2, function(x) { x * ((x - min(x)) / (max(x) - min(x))) })


  ##  Group SOMs

  if (length(unique(group.labels)) > 1) # mean group metagenes
  {
    group.metadata <<-
      do.call(cbind, by(t(metadata), group.labels, colMeans))[,unique(group.labels)]

    loglog.group.metadata <<-
      do.call(cbind, by(t(loglog.metadata), group.labels, colMeans))[,unique(group.labels)]

    WAD.group.metadata <<-
      do.call(cbind, by(t(WAD.metadata), group.labels, colMeans))[,unique(group.labels)]
  }


  ## set up SOM dependent variables
  genes.coordinates <<- apply(som.result$visual[,c(1,2)]+1, 1, paste, collapse=" x ")
  names(genes.coordinates) <<- rownames(indata)

  som.nodes <<- (som.result$visual[,"x"] + 1) + som.result$visual[,"y"] * preferences$dim.1stLvlSom
  names(som.nodes) <<- rownames(indata)

  return(T)
}
