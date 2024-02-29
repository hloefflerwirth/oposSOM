workspace.check <- function(env)
{
  if(missing("env")) stop("environment missing!")
  
  cat("Perform Workspace Check\n***********************\n"); flush.console()
  
  # Workspace objects and preferences list

  std <- oposSOM::opossom.new()
  std$preferences <- c(std$preferences,session.info="",started="",system.info="")
  
  
  if( length( setdiff( ls(env), ls(std) ) ) > 0 )
  {
    cat("spare objects in workspace:\n",paste( paste("\t",setdiff(ls(env),ls(std) )), collapse="\n" ),"\n"); flush.console()
  }
  
  if( length( setdiff( ls(std), ls(env) ) ) > 0 )
  {
    cat("missing objects in workspace:\n",paste( paste("\t",setdiff( ls(std), ls(env) )), collapse="\n" ),"\n"); flush.console()
  }
  
  if( length( setdiff( ls(env$preferences), ls(std$preferences) ) ) > 0 )
  {
    cat("spare parameters in workspace preferences:\n",paste( paste("\t",setdiff(ls(env$preferences),ls(std$preferences) )), collapse="\n" ),"\n"); flush.console()
  } 

  if( length( setdiff( ls(std$preferences), ls(env$preferences) ) ) > 0 )
  {
    cat("missing parameters in workspace preferences:\n",paste( paste("\t",setdiff( ls(std$preferences), ls(env$preferences) )), collapse="\n" ),"\n"); flush.console()
  }


  # Primary data check

  if (nrow(env$indata) == 0 || ncol(env$indata) == 0)
  {
    cat("indata: empty object\n"); flush.console()
  }
  if (any(is.na(rownames(env$indata))) || any(is.na(colnames(env$indata))))
  {
    cat("indata: NA in dimnames\n"); flush.console()
  }  
  if (any(rownames(env$indata)=="") || any(colnames(env$indata)==""))
  {
    cat("indata: \"\" in dimnames\n"); flush.console()
  }
  if (nrow(env$metadata) == 0 || ncol(env$metadata) == 0)
  {
    cat("metadata: empty object\n"); flush.console()
  }
  if (!all(colnames(env$metadata) == colnames(env$indata)))
  {
    cat("metadata: not the same samples as indata\n"); flush.console()
  }


  # som.result objects

  if (!"feature.BMU" %in% names(env$som.result) )
  {
    cat("som.result$feature.BMU: does not exist\n"); flush.console()
  }
  if (!all(names(env$som.result$feature.BMU) == rownames(env$indata)))
  {
    cat("som.result$feature.BMU: does not fit to features\n"); flush.console()
  }
  if (!"node.summary" %in% names(env$som.result) )
  {
    cat("som.result$node.summary: does not exist\n"); flush.console()
  } else
  if (colnames(env$som.result$node.summary)[3] != "n.features")
  {
    cat("som.result$node.summary: column n.features missing\n"); flush.console()
  } 
  
  
  # Data objects
  
  if (ncol(env$spot.list.overexpression$spotdata) != ncol(env$indata))
  {
    cat("spot.list.overexpression*: spotdata missing\n"); flush.console()
  }
  if (!"beta.statistic" %in% names(env$spot.list.overexpression$spots[[1]]))
  {
    cat("spot.list.*: beta scores missing\n"); flush.console()
  }  

  if( any( colnames(env$spot.list.correlation$spotdata)!=colnames(env$indata) ) )
  {
    cat("spot.list.correlation: spotdata columns do not fit indata columns\n"); flush.console()
  }
  if( any( colnames(env$spot.list.dmap$spotdata)!=colnames(env$indata) ) )
  {
    cat("spot.list.dmap: spotdata columns do not fit indata columns\n"); flush.console()
  }
  if( any( colnames(env$spot.list.group.overexpression$spotdata)!=colnames(env$indata) ) )
  {
    cat("spot.list.group.overexpression: spotdata columns do not fit indata columns\n"); flush.console()
  }
  if( any( colnames(env$spot.list.kmeans$spotdata)!=colnames(env$indata) ) )
  {
    cat("spot.list.kmeans: spotdata columns do not fit indata columns\n"); flush.console()
  }
  if( any( colnames(env$spot.list.overexpression$spotdata)!=colnames(env$indata) ) )
  {
    cat("spot.list.overexpression: spotdata columns do not fit indata columns\n"); flush.console()
  }
  if( any( colnames(env$spot.list.underexpression$spotdata)!=colnames(env$indata) ) )
  {
    cat("spot.list.underexpression: spotdata columns do not fit indata columns\n"); flush.console()
  }
  
  if( any( colnames(env$samples.GSZ.scores)!=colnames(env$indata) ) )
  {
    cat("samples.GSZ.scores: columns do not fit indata columns\n"); flush.console()
  } 
  

  # groups

  if (length(env$group.labels) != ncol(env$indata) || length(names(env$group.labels))!=ncol(env$indata) ||  !all(names(env$group.labels) == colnames(env$indata)))
  {
    cat("group.labels: does not fit to indata\n"); flush.console()
  }
  if (length(env$group.colors) != ncol(env$indata) || length(names(env$group.colors))!=ncol(env$indata) || !all(names(env$group.colors) == colnames(env$indata)))
  {
    cat("group.colors: does not fit to indata\n"); flush.console()
  }
  if (unique(substr(env$group.colors, 1, 1))[1] != "#")
  {
    cat("group.colors: not converted into #RGB format\n"); flush.console()
  }
  if (!all(env$group.colors %in% env$groupwise.group.colors))
  {
    cat("groupwise.group.colors: does not fit to group.colors\n"); flush.console()
  }
  if (!all(names(env$groupwise.group.colors) == unique(env$group.labels)))
  {
    cat("groupwise.group.colors: does not fit to group.labels\n"); flush.console()
  }
  if (unique(substr(env$groupwise.group.colors, 1, 1))[1] != "#")
  {
    cat("groupwise.group.colors: not converted into #RGB format\n"); flush.console()
  }
  if (!all(names(env$group.silhouette.coef) == colnames(env$indata)))
  {
    cat("group.silhouette.coef: does not fit to samples\n"); flush.console()
  }

   
  # Info objects
   
  std$gene.info <- list( ids = NULL, names = NULL, descriptions = NULL, coordinates = NULL,
                          chr.name = NULL, chr.band = NULL, chr.start = NULL, ensembl.mapping=NULL )
  
  for( x in names(std$gene.info) )
  { 
    if (! x %in% names(env$gene.info) )
    {
      cat("gene.info$",x,"not found in gene.info\n"); flush.console()
    } else
    if( x!="ensembl.mapping")
      if ( length(names(env$gene.info[[x]])) != nrow(env$indata) || !all(names(env$gene.info[[x]]) == rownames(env$indata)))
      {
        cat("gene.info$",x,": does not fit to features\n"); flush.console()
      }
  }
  

  # Geneset objects

  if (!is.null(env$gs.def.list))
  {
    unique.ens.ids <- unique(env$gene.info$ensembl.mapping$ensembl_gene_id)
    n.bad.sets = sum(!sapply(env$gs.def.list,function(x)all(x$Genes%in%unique.ens.ids)) )
    if (n.bad.sets>0)
    {
      cat("gs.def.list: genes of",n.bad.sets,"sets not in unique.protein.ids\n"); flush.console()
    }
  }


  # functions

  if (class(c) != "function")
  {
    cat("basic function: c overwritten\n"); flush.console()
  }
  if (class(t) != "function")
  {
    cat("basic function: t overwritten\n"); flush.console()
  }
  if ( F != FALSE )
  {
    cat("basic constant: F overwritten\n"); flush.console()
  }
  if ( T != TRUE )
  {
    cat("basic constant: T overwritten\n"); flush.console()
  }
  

  cat("\nThats all folks!\n****************\n"); flush.console()
}