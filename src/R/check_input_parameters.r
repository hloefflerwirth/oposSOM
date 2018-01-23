pipeline.checkInputParameters <- function()
{
  #### check preferences ####
  if (!is.character(preferences$dataset.name))
  {
    util.warn("Invalid value of \"dataset.name\". Using \"Unnamed\"")
    preferences$dataset.name <<- "Unnamed"
  }
  
  if ( preferences$dim.1stLvlSom!="auto" && !is.numeric(preferences$dim.1stLvlSom) || preferences$dim.1stLvlSom < 1)
  {
    util.warn("Invalid value of \"dim.1stLvlSom\". Using size recommendation")
    preferences$dim.1stLvlSom <<- "auto"
  }
  
  if (!is.numeric(preferences$dim.2ndLvlSom) || preferences$dim.2ndLvlSom < 1)
  {
    util.warn("Invalid value of \"dim.2ndLvlSom\". Using 20")
    preferences$dim.2ndLvlSom <<- 20
  }
  
  if (!is.numeric(preferences$training.extension) ||
      preferences$training.extension < 1 ||
      preferences$training.extension > 100)
  {
    util.warn("Invalid value of \"training.extension\". Using 1")
    preferences$training.extension <<- 1
  }
  
  if (!is.numeric(preferences$rotate.SOM.portraits) ||
      preferences$rotate.SOM.portraits < 0 ||
      preferences$rotate.SOM.portraits > 4)
  {
    util.warn("Invalid value of \"rotate.SOM.portraits\". Using 0")
    preferences$rotate.SOM.portraits <<- 0
  }
  
  if (!is.logical(preferences$flip.SOM.portraits))
  {
    util.warn("Invalid value of \"flip.SOM.portraits\". Using FALSE")
    preferences$flip.SOM.portraits <<- FALSE
  }
  
  if (!is.character(preferences$database.biomart))
  {
    util.warn("Invalid value of \"database.biomart\". Using \"\"")
    preferences$database.biomart <<- ""
  }
  
  if (!is.character(preferences$database.host))
  {
    util.warn("Invalid value of \"database.host\". Using \"\"")
    preferences$database.host <<- ""
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
  
  if (!is.list(preferences$activated.modules))
  {
    util.warn("Invalid value of \"activated.modules\". Using all analysis modules")
    preferences$activated.modules = list( "reporting" = TRUE,
                              "primary.analysis" = TRUE, 
                              "sample.similarity.analysis" = TRUE,
                              "geneset.analysis" = TRUE, 
                              "geneset.analysis.exact" = FALSE,
                              "group.analysis" = TRUE,
                              "difference.analysis" = TRUE )
  } else
  {
    if (!is.logical(preferences$activated.modules$reporting))
    {
      util.warn("Invalid value of \"activated.modules$reporting\". Using TRUE")
      preferences$activated.modules$reporting <<- TRUE
    }
    if (!is.logical(preferences$activated.modules$primary.analysis))
    {
      util.warn("Invalid value of \"activated.modules$primary.analysis\". Using TRUE")
      preferences$activated.modules$primary.analysis <<- TRUE
    } else
    if (!preferences$activated.modules$primary.analysis && is.null(env$som.result))
    {
      util.warn("No primary analysis perfomed yet. Setting \"activated.modules$primary.analysis\" to TRUE")
      preferences$activated.modules$primary.analysis <<- TRUE
    }
    if (!is.logical(preferences$activated.modules$sample.similarity.analysis))
    {
      util.warn("Invalid value of \"activated.modules$sample.similarity.analysis\". Using TRUE")
      preferences$activated.modules$sample.similarity.analysis <<- TRUE
    }
    if (!is.logical(preferences$activated.modules$geneset.analysis))
    {
      util.warn("Invalid value of \"activated.modules$geneset.analysis\". Using TRUE")
      preferences$activated.modules$geneset.analysis <<- TRUE
    }
    if (!is.logical(preferences$activated.modules$geneset.analysis.exact))
    {
      util.warn("Invalid value of \"activated.modules$geneset.analysis.exact\". Using FALSE")
      preferences$activated.modules$geneset.analysis.exact <<- FALSE
    }    
    if (!is.logical(preferences$activated.modules$group.analysis))
    {
      util.warn("Invalid value of \"activated.modules$group.analysis\". Using TRUE")
      preferences$activated.modules$group.analysis <<- TRUE
    }
    if (!is.logical(preferences$activated.modules$difference.analysis))
    {
      util.warn("Invalid value of \"activated.modules$difference.analysis\". Using TRUE")
      preferences$activated.modules$difference.analysis <<- TRUE
    }
  }
  
  if (!is.character(preferences$standard.spot.modules) || length(preferences$standard.spot.modules)!=1 ||
      !preferences$standard.spot.modules %in% c("overexpression","underexpression","kmeans","correlation","group.overexpression","dmap") )
  {
    util.warn("Invalid value of \"standard.spot.modules\". Using \"dmap\"")
    preferences$standard.spot.modules <<- "dmap"
  }
  
  if (!is.numeric(preferences$spot.coresize.modules) ||
      preferences$spot.coresize.modules < 1 ||
      preferences$spot.coresize.modules > 20)
  {
    util.warn("Invalid value of \"spot.coresize.modules\". Using 3")
    preferences$spot.coresize.modules <<- 3
  }
  
  if (!is.numeric(preferences$spot.threshold.modules) ||
      preferences$spot.threshold.modules <= 0 ||
      preferences$spot.threshold.modules >= 1)
  {
    util.warn("Invalid value of \"spot.threshold.modules\". Using 0.95")
    preferences$spot.threshold.modules <<- 0.95
  }
  
  if (!is.numeric(preferences$spot.coresize.groupmap) ||
      preferences$spot.coresize.groupmap < 1 ||
      preferences$spot.coresize.groupmap > 20)
  {
    util.warn("Invalid value of \"spot.coresize.groupmap\". Using 5")
    preferences$spot.coresize.groupmap <<- 5
  }
  
  if (!is.numeric(preferences$spot.threshold.groupmap) ||
      preferences$spot.threshold.groupmap <= 0 ||
      preferences$spot.threshold.groupmap >= 1)
  {
    util.warn("Invalid value of \"spot.threshold.groupmap\". Using 0.75")
    preferences$spot.threshold.groupmap <<- 0.75
  }
  
  if (!is.logical(preferences$feature.centralization))
  {
    util.warn("Invalid value of \"feature.centralization\". Using TRUE")
    preferences$feature.centralization <<- TRUE
  }
  
  if (!is.logical(preferences$sample.quantile.normalization))
  {
    util.warn("Invalid value of \"sample.quantile.normalization\". Using TRUE")
    preferences$sample.quantile.normalization <<- TRUE
  }
  
  if (!is.null(preferences$pairwise.comparison.list))
  {
    # translate sample names into indexes
    preferences$pairwise.comparison.list <<- lapply(preferences$pairwise.comparison.list,function(x)
    {
      if( class(x[[1]]) == "character" ) x[[1]] <- match(x[[1]],colnames(indata))
      if( class(x[[2]]) == "character" ) x[[2]] <- match(x[[2]],colnames(indata))
      return(x)
    })  
    
    # seek for corrupt sets
    empty.sets <- which( sapply(preferences$pairwise.comparison.list, function(x) min(sapply(x,function(y) length(na.omit(y))))) == 0 )
    if(length(empty.sets)>0)
    {   
      util.warn("Empty sample set found and removed from \"pairwise.comparison.list\".")
      preferences$pairwise.comparison.list <<- preferences$pairwise.comparison.list[-empty.sets] 
    }
    
    if (is.null(names(preferences$pairwise.comparison.list)) )
    {
      names(preferences$pairwise.comparison.list) <<- sapply(preferences$pairwise.comparison.list,function(x) paste(names(x), collapse=" vs ") )
    }
    names(preferences$pairwise.comparison.list)[which(names(preferences$pairwise.comparison.list)=="")] <<- 1:sum(names(preferences$pairwise.comparison.list)=="")   
  }
  
  
  #### check input data ####
  if (is.null(indata))
  {
    util.fatal("No indata supplied!")
    return(FALSE)
  }
  
  if (class(indata) == "ExpressionSet")
  {
    group.labels <<- as.character(pData(indata)$group.labels)
    group.colors <<- as.character(pData(indata)$group.colors)
    indata <<- assayData(indata)$exprs
  }
  
  if (class(indata) != "matrix" && (is.null(dim(indata)) || dim(indata) < 1))
  {
    util.fatal("Invalid indata! Provide a two-dimensional numerical matrix.")
    return(FALSE)
  }
  
  if (class(indata) != "matrix" ||
      mode(indata) != "numeric" )
      # storage.mode(indata) != "numeric")
  {
    rn <- rownames(indata)
    num.mode <- sapply(seq(ncol(indata)), function(i){ all(grepl("^-?[0-9\\.]+$", indata[,i])) })
    
    if( num.mode[1]==FALSE && all(num.mode[-1]==TRUE) ) # check if IDs are contained as first row
    {
      rn <- indata[,1]
      indata <<- indata[,-1]
      num.mode <- num.mode[-1]
      util.warn("Gene IDs adopted from first data column.")    
    }
    
    if( any(num.mode!=TRUE) ) # check if all columns contain numbers or convertable characters
    {
      util.fatal("Invalid indata! Provide a two-dimensional numerical matrix.")
      return(FALSE)
      
    } else
    {
      indata <<- apply(indata, 2, function(x){ as.numeric(as.vector(x)) })
      rownames(indata) <<- rn
      storage.mode(indata) <<- "numeric"
      util.warn("Indata converted to two-dimensional numerical matrix.")    
    }
  }
  
  if( length(group.labels)==1 && group.labels=="auto" )
  {
    group.labels <<- rep("auto",ncol(indata)) 
    names(group.labels) <<- colnames(indata)
  }
  
  const.cols <- which(apply(indata, 2, function(col) { diff(range(col)) == 0 }))
  
  if (length(const.cols) > 0)
  {
    indata <<- indata[,-const.cols]
    group.labels <<- group.labels[-const.cols]
    group.colors <<- group.colors[-const.cols]
    util.warn("Removed",length(const.cols),"constant columns from data set.")
  }
  
  const.rows <- which(apply(indata, 1, function(row) { diff(range(row)) == 0 }))
  
  if (length(const.rows) > 0)
  {
    indata <<- indata[-const.rows,]
    util.warn("Removed",length(const.rows),"constant rows from data set.")
  }
  
  if (length(rownames(indata)) == 0)
  {
    rownames(indata) <<- as.character(1:nrow(indata))
    preferences$activated.modules$geneset.analysis <<- FALSE
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
  if( "" %in% rownames(indata) )
  {
    indata <<- indata[which(rownames(indata)!=""),]
    util.warn("Removed genes with \"\" id from data")
  }
  
  na.rows <- which( apply(indata, 1, function(x) sum( is.na(x) | is.infinite(x) ) ) > 0 )
  
  if (length(na.rows) > 0)
  {
    indata <<- indata[-na.rows,]
    util.warn("Removed NAs or infinite values from data set")
  }
  
  if (preferences$dim.1stLvlSom == "auto")
  {
    n.sample.interval <- cut( ncol(indata), breaks=c(0,100,500,1000,5000,Inf), labels=c(1:5) )
    n.feature.interval <- cut( nrow(indata), breaks=c(0,1000,10000,Inf), labels=c(1:3) )
    recommendation <- matrix(c(seq(20,40,5),seq(30,50,5),seq(40,60,5)),nrow=3,byrow=TRUE)
    
    preferences$dim.1stLvlSom <<- recommendation[n.feature.interval,n.sample.interval]
    util.info("Recommended SOM size will be used:",preferences$dim.1stLvlSom,"x",preferences$dim.1stLvlSom) 
  }
  
  # check group.labels and group.colors
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
    for (sample in unique(colnames(indata)))
    {
      if (length(unique(group.labels[which(colnames(indata) == sample)])) > 1)
      {
        util.warn("Sample is in multiple groups:", sample)
        group.labels <<- NULL
        group.colors <<- NULL
        break
      }
    }
  }
  
  if (!is.null(group.labels))
  {
    group.labels <<- as.character(group.labels)
    names(group.labels) <<- colnames(indata)
    
    if (is.null(group.colors))
    {
      group.colors <<- rep("", ncol(indata))
      
      for (i in seq_along(unique(group.labels)))
      {
        group.colors[which(group.labels == unique(group.labels)[i])] <<-
          colorRampPalette(c("blue3", "blue", "lightblue", "green2", "gold", "red", "red3"))(length(unique(group.labels)))[i]
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
    group.labels <<- rep("auto",ncol(indata))
    names(group.labels) <<- colnames(indata)
    
    group.colors <<- rep("#000000", ncol(indata))
    names(group.colors) <<- colnames(indata)
  }
  
  groupwise.group.colors <<- group.colors[match(unique(group.labels), group.labels)]
  names(groupwise.group.colors) <<- unique(group.labels)
  
  
  # set color schemes
  if (!is.null(color.palette.portraits)) # check if given color palette is a valid function
  {
    if( length(environment(color.palette.portraits))!=3 || !all( c("colors","ramp") %in% ls(environment(color.palette.portraits)) ) )
    {
      util.warn("Invalid value of \"color.palette.portraits\". Using standard scheme")
      color.palette.portraits <<- colorRampPalette(c("darkblue","blue","lightblue","green2","yellow","red","darkred"))
    }
  } else
  {
    color.palette.portraits <<- colorRampPalette(c("darkblue","blue","lightblue","green2","yellow","red","darkred"))
  }
  
  if (!is.null(color.palette.heatmaps)) # check if given color palette is a valid function
  {
    if( length(environment(color.palette.heatmaps))!=3 || !all( c("colors","ramp") %in% ls(environment(color.palette.heatmaps)) ) )
    {
      util.warn("Invalid value of \"color.palette.heatmaps\". Using standard scheme")
      color.palette.heatmaps <<- colorRampPalette(c("blue4","blue","gray90","orange","red4"))
    }
  } else
  {
    color.palette.heatmaps <<- colorRampPalette(c("blue4","blue","gray90","orange","red4"))
  }
  
  if(preferences$activated.modules$primary.analysis)
  {
    files.name <<- preferences$dataset.name
    while (file.exists(paste(files.name, ".RData", sep=""))) {
      files.name <<- paste(files.name, "+", sep="")
    }
    
    output.paths <<-
      c("CSV"=paste(files.name, "- Results/CSV Sheets"),
        "Summary Sheets Samples"=paste(files.name, "- Results/Summary Sheets - Samples") )
  } 
  
  return(TRUE)
}