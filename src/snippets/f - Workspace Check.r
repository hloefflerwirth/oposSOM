
  attach(env)

  cat("Perform Workspace Check\n\n"); flush.console()
  
  
    # standard preferences
        
    if (! ("geneset.analysis.samplespots" %in% names(preferences)))
    {
      cat("preferences: geneset.analysis.samplespots not set correctly\n"); flush.console()
    } 
  

    # primary data
  
    if (nrow(indata) == 0 || ncol(indata) == 0)
    {
      cat("indata: empty object\n"); flush.console()
    }
    if (any(is.na(rownames(indata))) || any(is.na(colnames(indata))))  
    {
      cat("indata: NA in dimnames\n"); flush.console()
    }      
    if (nrow(metadata) == 0 || ncol(metadata) == 0)
    {
      cat("metadata: empty object\n"); flush.console()
    }  
    if (!all(colnames(metadata) == colnames(indata)))
    {
      cat("metadata: not the same samples as indata\n"); flush.console()
    }  
  
  
    # secondary objects
  
    if (!all(names(spot.list.samples) == colnames(indata)))
    {      
      cat("GS.infos.samples: not the same samples as indata\n"); flush.console()        
    }  
    if (!all(names(som.nodes) == rownames(indata)))
    {
      cat("som.nodes: does not fit to features\n"); flush.console()            
    }  
  
  
  
  
  
    # groups
  
    if (length(group.labels) != ncol(indata) || length(names(group.labels))!=ncol(indata) ||  !all(names(group.labels) == colnames(indata)))
    {
      cat("group.labels: does not fit to indata\n"); flush.console()
    }  
    if (length(group.colors) != ncol(indata) || length(names(group.colors))!=ncol(indata) || !all(names(group.colors) == colnames(indata)))
    {
      cat("group.colors: does not fit to indata\n"); flush.console()
    }  
    if (unique(substr(group.colors, 1, 1))[1] != "#")
    {
      cat("group.colors: not converted into #RGB format\n"); flush.console()
    }  
    if (!all(group.colors %in% groupwise.group.colors))
    {
      cat("groupwise.group.colors: does not fit to group.colors\n"); flush.console()            
    }
    if (!all(names(groupwise.group.colors) == unique(group.labels)))
    {
      cat("groupwise.group.colors: does not fit to group.labels\n"); flush.console()            
    }
    if (unique(substr(groupwise.group.colors, 1, 1))[1] != "#")
    {
      cat("groupwise.group.colors: not converted into #RGB format\n"); flush.console()
    }
  
    if (!all(names(group.silhouette.coef) == colnames(indata)))
    {
      cat("group.bootstrap.score: does not fit to samples\n"); flush.console()            
    }  
  
    
  
  
    # statistic objects
  
    if (!all(colnames(WAD.g.m) == colnames(indata)) || !all(rownames(WAD.g.m) == rownames(indata)))
    {
      cat("WAD.g.m: does not fit to indata\n"); flush.console()            
    }      
    if (!all(colnames(t.g.m) == colnames(indata)) || !all(rownames(t.g.m) == rownames(indata)))
    {
      cat("t.g.m: does not fit to indata\n"); flush.console()            
    }  
    if (!all(colnames(p.g.m) == colnames(indata)) || !all(rownames(p.g.m) == rownames(indata)))
    {
      cat("p.g.m: does not fit to indata\n"); flush.console()            
    }
    if (!all(colnames(fdr.g.m) == colnames(indata)) || !all(rownames(fdr.g.m) == rownames(indata)))
    {
      cat("fdr.g.m: does not fit to indata\n"); flush.console()            
    }
    if (!all(colnames(Fdr.g.m) == colnames(indata)) || !all(rownames(Fdr.g.m) == rownames(indata)))
    {
      cat("Fdr.g.m: does not fit to indata\n"); flush.console()            
    }

  
    # info objects
  
    if (!all(names(gene.descriptions) == rownames(indata)))
    {
      cat("gene.descriptions: does not fit to features\n"); flush.console()            
    }      
    if (!all(names(gene.names) == rownames(indata)))
    {
      cat("gene.names: does not fit to features\n"); flush.console()            
    }
    if (!all(names(gene.ids) %in% rownames(indata)))
    {
      cat("gene.ids: does not fit to features\n"); flush.console()            
    }  
    if (!all(names(gene.positions) %in% rownames(indata)))
    {
      cat("gene.positions: does not fit to features\n"); flush.console()            
    }    
    if (!all(names(gene.coordinates) == rownames(indata)))
    {
      cat("gene.coordinates: does not fit to features\n"); flush.console()            
    }  
    
  
  
    # Pipeline data structures
  
    if (ncol(spot.list.overexpression$spotdata) != ncol(indata))
    {
      cat("spot.list.overexpression*: spotdata missing\n"); flush.console()            
    }  
    if (!"beta.statistic" %in% names(spot.list.overexpression$spots[[1]]))
    {
      cat("spot.list.*: beta scores missing\n"); flush.console()            
    }  
    if ( !exists("spot.list.dmap") )
    {
      cat("spot.list.dmap missing\n"); flush.console()            
    }  
  
  
    # basic functions sometimes overwritten
  
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
  

  cat("\nThats all folks!\n\n"); flush.console()
  
