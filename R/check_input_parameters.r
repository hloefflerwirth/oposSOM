pipeline.checkInputParameters <- function(env)
{
  #### check preferences ####
  if (!is.character(env$preferences$dataset.name))
  {
    util.warn("Invalid value of \"dataset.name\". Using \"Unnamed\"")
    env$preferences$dataset.name <- "Unnamed"
  }
  
  if ( env$preferences$dim.1stLvlSom!="auto" && !is.numeric(env$preferences$dim.1stLvlSom) || env$preferences$dim.1stLvlSom < 1)
  {
    util.warn("Invalid value of \"dim.1stLvlSom\". Using size recommendation")
    env$preferences$dim.1stLvlSom <- "auto"
  }
  
  if (!is.numeric(env$preferences$dim.2ndLvlSom) || env$preferences$dim.2ndLvlSom < 1)
  {
    util.warn("Invalid value of \"dim.2ndLvlSom\". Using 20")
    env$preferences$dim.2ndLvlSom <- 20
  }
  
  if (!is.numeric(env$preferences$training.extension) ||
      env$preferences$training.extension < 1 ||
      env$preferences$training.extension > 100)
  {
    util.warn("Invalid value of \"training.extension\". Using 1")
    env$preferences$training.extension <- 1
  }
  
  if (!is.numeric(env$preferences$rotate.SOM.portraits) ||
      env$preferences$rotate.SOM.portraits < 0 ||
      env$preferences$rotate.SOM.portraits > 4)
  {
    util.warn("Invalid value of \"rotate.SOM.portraits\". Using 0")
    env$preferences$rotate.SOM.portraits <- 0
  }
  
  if (!is.logical(env$preferences$flip.SOM.portraits))
  {
    util.warn("Invalid value of \"flip.SOM.portraits\". Using FALSE")
    env$preferences$flip.SOM.portraits <- FALSE
  }
  
  if (!is.character(env$preferences$database.biomart))
  {
    util.warn("Invalid value of \"database.biomart\". Using \"\"")
    env$preferences$database.biomart <- ""
  }
  
  if (!is.character(env$preferences$database.host))
  {
    util.warn("Invalid value of \"database.host\". Using \"\"")
    env$preferences$database.host <- ""
  }
  
  if (!is.character(env$preferences$database.dataset))
  {
    util.warn("Invalid value of \"database.dataset\". Using \"\"")
    env$preferences$database.dataset <- ""
  }
  
  if (!is.character(env$preferences$database.id.type))
  {
    util.warn("Invalid value of \"database.id.type\". Using \"\"")
    env$preferences$database.id.type <- ""
  }
  
  if (!is.list(env$preferences$activated.modules))
  {
    util.warn("Invalid value of \"activated.modules\". Using all analysis modules")
    env$preferences$activated.modules = list( "reporting" = TRUE,
                              "primary.analysis" = TRUE, 
                              "sample.similarity.analysis" = TRUE,
                              "geneset.analysis" = TRUE, 
                              "psf.analysis" = TRUE,
                              "group.analysis" = TRUE,
                              "difference.analysis" = TRUE )
  } else
  {
    if (!is.logical(env$preferences$activated.modules$reporting))
    {
      util.warn("Invalid value of \"activated.modules$reporting\". Using TRUE")
      env$preferences$activated.modules$reporting <- TRUE
    }
    if (!is.logical(env$preferences$activated.modules$primary.analysis))
    {
      util.warn("Invalid value of \"activated.modules$primary.analysis\". Using TRUE")
      env$preferences$activated.modules$primary.analysis <- TRUE
    } else
    if (!env$preferences$activated.modules$primary.analysis && is.null(env$som.result))
    {
      util.warn("No primary analysis perfomed yet. Setting \"activated.modules$primary.analysis\" to TRUE")
      env$preferences$activated.modules$primary.analysis <- TRUE
    }
    if (!is.logical(env$preferences$activated.modules$sample.similarity.analysis))
    {
      util.warn("Invalid value of \"activated.modules$sample.similarity.analysis\". Using TRUE")
      env$preferences$activated.modules$sample.similarity.analysis <- TRUE
    }
    if (!is.logical(env$preferences$activated.modules$geneset.analysis))
    {
      util.warn("Invalid value of \"activated.modules$geneset.analysis\". Using TRUE")
      env$preferences$activated.modules$geneset.analysis <- TRUE
    }
    if (!is.logical(env$preferences$activated.modules$psf.analysis))
    {
      util.warn("Invalid value of \"activated.modules$psf.analysis\". Using TRUE")
      env$preferences$activated.modules$psf.analysis <- TRUE
    }    
    if (!is.logical(env$preferences$activated.modules$group.analysis))
    {
      util.warn("Invalid value of \"activated.modules$group.analysis\". Using TRUE")
      env$preferences$activated.modules$group.analysis <- TRUE
    }
    if (!is.logical(env$preferences$activated.modules$difference.analysis))
    {
      util.warn("Invalid value of \"activated.modules$difference.analysis\". Using TRUE")
      env$preferences$activated.modules$difference.analysis <- TRUE
    }
  }
  
  if (!is.character(env$preferences$standard.spot.modules) || length(env$preferences$standard.spot.modules)!=1 ||
      !env$preferences$standard.spot.modules %in% c("overexpression","underexpression","kmeans","correlation","group.overexpression","dmap") )
  {
    util.warn("Invalid value of \"standard.spot.modules\". Using \"dmap\"")
    env$preferences$standard.spot.modules <- "dmap"
  }
  
  if (!is.numeric(env$preferences$spot.coresize.modules) ||
      env$preferences$spot.coresize.modules < 1 ||
      env$preferences$spot.coresize.modules > 20)
  {
    util.warn("Invalid value of \"spot.coresize.modules\". Using 3")
    env$preferences$spot.coresize.modules <- 3
  }
  
  if (!is.numeric(env$preferences$spot.threshold.modules) ||
      env$preferences$spot.threshold.modules <= 0 ||
      env$preferences$spot.threshold.modules >= 1)
  {
    util.warn("Invalid value of \"spot.threshold.modules\". Using 0.95")
    env$preferences$spot.threshold.modules <- 0.95
  }
  
  if (!is.numeric(env$preferences$spot.coresize.groupmap) ||
      env$preferences$spot.coresize.groupmap < 1 ||
      env$preferences$spot.coresize.groupmap > 20)
  {
    util.warn("Invalid value of \"spot.coresize.groupmap\". Using 5")
    env$preferences$spot.coresize.groupmap <- 5
  }
  
  if (!is.numeric(env$preferences$spot.threshold.groupmap) ||
      env$preferences$spot.threshold.groupmap <= 0 ||
      env$preferences$spot.threshold.groupmap >= 1)
  {
    util.warn("Invalid value of \"spot.threshold.groupmap\". Using 0.75")
    env$preferences$spot.threshold.groupmap <- 0.75
  }
  
  if (!is.logical(env$preferences$feature.centralization))
  {
    util.warn("Invalid value of \"feature.centralization\". Using TRUE")
    env$preferences$feature.centralization <- TRUE
  }
  
  if (!is.logical(env$preferences$sample.quantile.normalization))
  {
    util.warn("Invalid value of \"sample.quantile.normalization\". Using TRUE")
    env$preferences$sample.quantile.normalization <- TRUE
  }
  
  if (!is.null(env$preferences$pairwise.comparison.list))
  {
    # translate sample names into indexes
    env$preferences$pairwise.comparison.list <- lapply(env$preferences$pairwise.comparison.list,function(x)
    {
      if( is(x[[1]],"character") ) x[[1]] <- match(x[[1]],colnames(env$indata))
      if( is(x[[2]],"character") ) x[[2]] <- match(x[[2]],colnames(env$indata))
      return(x)
    })  
    
    # seek for corrupt sets
    empty.sets <- which( sapply(env$preferences$pairwise.comparison.list, function(x) min(sapply(x,function(y) length(na.omit(y))))) == 0 )
    if(length(empty.sets)>0)
    {   
      util.warn("Empty sample set found and removed from \"pairwise.comparison.list\".")
      env$preferences$pairwise.comparison.list <- preferences$pairwise.comparison.list[-empty.sets] 
    }
    
    if (is.null(names(env$preferences$pairwise.comparison.list)) )
    {
      names(env$preferences$pairwise.comparison.list) <- sapply(env$preferences$pairwise.comparison.list,function(x) paste(names(x), collapse=" vs ") )
    }
    names(env$preferences$pairwise.comparison.list)[which(names(env$preferences$pairwise.comparison.list)=="")] <- 1:sum(names(env$preferences$pairwise.comparison.list)=="")   
  }
  
  if (!is.numeric(env$preferences$max.cores) ||
      env$preferences$max.cores < 1 ||
      env$preferences$max.cores > detectCores()-1 )
  {
    env$preferences$max.cores <- detectCores()-1
    util.warn( paste( "Invalid value of \"max.cores\". Using", detectCores()-1 ) )
  }
  
  #### check input data ####
  if (is.null(env$indata))
  {
    util.fatal("No indata supplied!")
    env$passedInputChecking <- FALSE
    return(env)
  }
  
  if (is(env$indata,"ExpressionSet"))
  {
    env$group.labels <- as.character(pData(env$indata)$group.labels)
    env$group.colors <- as.character(pData(env$indata)$group.colors)
    env$indata <- assayData(env$indata)$exprs
  }
  
  if (!is(env$indata,"matrix") && (is.null(dim(env$indata)) || dim(env$indata) < 1))
  {
    util.fatal("Invalid indata! Provide a two-dimensional numerical matrix.")
    env$passedInputChecking <- FALSE
    return(env)
  }
  
  if (!is(env$indata,"matrix") ||
      mode(env$indata) != "numeric" )
      # storage.mode(indata) != "numeric")
  {
    rn <- rownames(env$indata)
    num.mode <- sapply(seq(ncol(env$indata)), function(i){ all(grepl("^-?[0-9\\.]+$", env$indata[,i])) })
    
    if( num.mode[1]==FALSE && all(num.mode[-1]==TRUE) ) # check if IDs are contained as first row
    {
      rn <- env$indata[,1]
      env$indata <- env$indata[,-1]
      num.mode <- num.mode[-1]
      util.warn("Gene IDs adopted from first data column.")    
    }
    
    if( any(num.mode!=TRUE) ) # check if all columns contain numbers or convertable characters
    {
      util.fatal("Invalid indata! Provide a two-dimensional numerical matrix.")
      env$passedInputChecking <- FALSE
      return(env)
      
    } else
    {
      env$indata <- apply(env$indata, 2, function(x){ as.numeric(as.vector(x)) })
      rownames(env$indata) <- rn
      storage.mode(env$indata) <- "numeric"
      util.warn("Indata converted to two-dimensional numerical matrix.")    
    }
  }
  
  if( length(env$group.labels)==1 && env$group.labels=="auto" )
  {
    env$group.labels <- rep("auto",ncol(env$indata)) 
    names(env$group.labels) <- colnames(env$indata)
  }
  
  const.cols <- which(apply(env$indata, 2, function(col) { diff(range(col)) == 0 }))
  
  if (length(const.cols) > 0)
  {
    env$indata <- env$indata[,-const.cols]
    env$group.labels <- env$group.labels[-const.cols]
    env$group.colors <- env$group.colors[-const.cols]
    util.warn("Removed",length(const.cols),"constant columns from data set.")
  }
  
  const.rows <- which(apply(env$indata, 1, function(row) { diff(range(row)) == 0 }))
  
  if (length(const.rows) > 0)
  {
    env$indata <- env$indata[-const.rows,]
    util.warn("Removed",length(const.rows),"constant rows from data set.")
  }
  
  if (length(rownames(env$indata)) == 0)
  {
    rownames(env$indata) <- as.character(1:nrow(env$indata))
    env$preferences$activated.modules$geneset.analysis <- FALSE
    util.warn("No rownames found. Set them to 1,2,3,4...")
  }
  
  if (length(colnames(env$indata)) == 0)
  {
    colnames(env$indata) <- paste("Sample", c(1:ncol(env$indata)))
    util.warn("No colnames found. Set them to 1,2,3,4...")
  }
  
  if (any(duplicated(rownames(env$indata))))
  {
    env$indata <- do.call(rbind, by(env$indata, rownames(env$indata), colMeans))[unique(rownames(env$indata)),]
    util.warn("Duplicate rownames. Averaged multiple features")
  }
  if( "" %in% rownames(env$indata) )
  {
    env$indata <- env$indata[which(rownames(env$indata)!=""),]
    util.warn("Removed genes with \"\" id from data")
  }
  
  na.rows <- which( apply(env$indata, 1, function(x) sum( is.na(x) | is.infinite(x) ) ) > 0 )
  
  if (length(na.rows) > 0)
  {
    env$indata <- env$indata[-na.rows,]
    util.warn("Removed NAs or infinite values from data set")
  }
  
  if (env$preferences$dim.1stLvlSom == "auto")
  {
    n.sample.interval <- cut( ncol(env$indata), breaks=c(0,100,500,1000,5000,Inf), labels=c(1:5) )
    n.feature.interval <- cut( nrow(env$indata), breaks=c(0,1000,10000,Inf), labels=c(1:3) )
    recommendation <- matrix(c(seq(20,40,5),seq(30,50,5),seq(40,60,5)),nrow=3,byrow=TRUE)
    
    env$preferences$dim.1stLvlSom <- recommendation[n.feature.interval,n.sample.interval]
    util.info("Recommended SOM size will be used:",env$preferences$dim.1stLvlSom,"x",env$preferences$dim.1stLvlSom) 
  }
  
  # check group.labels and group.colors
  if ((!is.null(env$group.labels) && length(env$group.labels) != ncol(env$indata)) ||
      (!is.null(env$group.colors) && length(env$group.colors) != ncol(env$indata)))
  {
    env$group.labels <- NULL
    env$group.colors <- NULL
    util.warn("Group assignment doesnt fit number of samples")
  }
  
  if (!is.null(env$group.labels) && max(table(env$group.labels)) == 1)
  {
    env$group.labels <- NULL
    env$group.colors <- NULL
    util.warn("Each sample has an own group")
  }
  
  if (!is.null(env$group.labels))
  {
    for (sample in unique(colnames(env$indata)))
    {
      if (length(unique(env$group.labels[which(colnames(env$indata) == sample)])) > 1)
      {
        util.warn("Sample is in multiple groups:", sample)
        env$group.labels <- NULL
        env$group.colors <- NULL
        break
      }
    }
  }
  
  if (!is.null(env$group.labels))
  {
    env$group.labels <- as.character(env$group.labels)
    names(env$group.labels) <- colnames(env$indata)
    
    if (is.null(env$group.colors))
    {
      env$group.colors <- rep("", ncol(env$indata))
      
      for (i in seq_along(unique(env$group.labels)))
      {
        env$group.colors[which(env$group.labels == unique(env$group.labels)[i])] <- color.palette.discrete(length(unique(env$group.labels)))[i]
      }
    }
    
    # catch userdefined group.colors --> convert to #hex
    if (length(unique(substr(env$group.colors, 1, 1)) > 1) || unique(substr(env$group.colors, 1, 1))[1] != "#")
    {
      env$group.colors <- apply(col2rgb(env$group.colors), 2, function(x) { rgb(x[1]/255, x[2]/255, x[3]/255) })
    }
    names(env$group.colors) <- colnames(env$indata)
  } else
  {
    env$group.labels <- rep("auto",ncol(env$indata))
    names(env$group.labels) <- colnames(env$indata)
    
    env$group.colors <- rep("#000000", ncol(env$indata))
    names(env$group.colors) <- colnames(env$indata)
  }
  
  env$groupwise.group.colors <- env$group.colors[match(unique(env$group.labels), env$group.labels)]
  names(env$groupwise.group.colors) <- unique(env$group.labels)
  
  
  # set color schemes
  if (!is.null(env$color.palette.portraits)) # check if given color palette is a valid function
  {
    if( length(environment(env$color.palette.portraits))!=3 || !all( c("colors","ramp") %in% ls(environment(env$color.palette.portraits)) ) )
    {
      util.warn("Invalid value of \"color.palette.portraits\". Using standard scheme")
      env$color.palette.portraits <- if(env$preferences$colorblindsave.portraits) color.palette.portraits.cbs else color.palette.portraits
    }
  } else
  {
    env$color.palette.portraits <- if(env$preferences$colorblindsave.portraits) color.palette.portraits.cbs else color.palette.portraits
  }
  
  if (!is.null(env$color.palette.heatmaps)) # check if given color palette is a valid function
  {
    if( length(environment(env$color.palette.heatmaps))!=3 || !all( c("colors","ramp") %in% ls(environment(env$color.palette.heatmaps)) ) )
    {
      util.warn("Invalid value of \"color.palette.heatmaps\". Using standard scheme")
      env$color.palette.heatmaps <- color.palette.heatmaps
    }
  } else
  {
    env$color.palette.heatmaps <- color.palette.heatmaps
  }
  
  if(env$preferences$activated.modules$primary.analysis)
  {
    env$files.name <- env$preferences$dataset.name
    while (file.exists(paste(env$files.name, "- Results"))) {
      env$files.name <- paste(env$files.name, "+", sep="")
    }
    
    env$output.paths <-
      c("CSV"="CSV Sheets","Summary Sheets Samples"="Summary Sheets - Samples")
    
    if( !grepl("german",sessionInfo()$locale,ignore.case=TRUE) &&
        !grepl("bioinf.uni-leipzig.de",Sys.info()["nodename"],ignore.case=TRUE) )
      env$csv.function <- write.csv
  } 
  
  if( is.null(env$preferences$activated.modules$largedata.mode) )
  {
	env$preferences$activated.modules$largedata.mode <- ( ncol(env$indata) >= 1000 )
	
	if( ncol(env$indata) >= 1000 )
	{
		util.warn("Data has >1000 samples. 'Large data mode' activated!")
	}     
  }
  
  util.info("Data dimension:",nrow(env$indata),"features in ",ncol(env$indata),"samples")
  env$passedInputChecking <- TRUE
  return(env)
}