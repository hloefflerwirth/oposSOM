


pathway.expression.mapping <- function( g, fc.m )
{
  g@nodeData@data <- lapply(g@nodeData@data, function(x, fc.m)
  {
    if( sum(x$genes%in%names(fc.m)) > 0)
    {
      x$expression <- mean(fc.m[x$genes],na.rm=T)
    }
    else
    {
      x$expression <- 1
    }
    return(x)
  }, fc.m )

  return(g)
}




psf.flow <- function (g, node.ordering, sink.nodes, nodeIDs, split=TRUE)
{
  node.order <- node.ordering$node.order
  node.rank <- node.ordering$node.rank
  nods <- names(node.order)
  symb.exprs <- vector("list")
  eval.exprs <- vector("list")
  E <- data.frame(sapply(g@nodeData@data[nods],function(x)x$expression),
                 row.names=nods)
  I <- matrix(data=NA, nrow=length(nods), ncol=length(nods), dimnames=list(nods,nods) )

  for (node in nods)
  {
    l <- length(edgeL(g)[[node]]$edges)
    for (e in 1:l) {
      from <- node
      to <- nodes(g)[edgeL(g)[node][[1]]$edges][e]
      if (!is.na(to) && length(to) > 0) {
        impact <- edgeData(g, from = from, to = to,
                                 attr = "impact")[[1]]
        I[from, to] <- impact
      }
    }
  }

  for (i in 1:length(nods))
  {
    node <- nods[i]
    parent.nodes <- inEdges(node, g)[[1]]

    if (length(parent.nodes) > 0)
    {
      node.exp <- E[i, 1]
      in.signal <- unlist(nodeData(g, parent.nodes,
                                          attr = "signal"))
      pi <- which(nods %in% parent.nodes)
      node.signal <- nodeData(g, nods[i], attr = "signal")[[1]]
      impact <- I[parent.nodes, nods[i]]
      if(split){
        proportion = in.signal/sum(in.signal)
        node.signal <- sum((proportion*node.exp)*(in.signal^impact))
      } else {
        #Returns the product of signals without splitting - is for updateing, but applies only in this special case where all the rules are s*t, or 1/s*t
        proportion = 1
        node.signal <- prod((proportion*node.exp)*(in.signal^impact))
      }
      nodeData(g, nods[i], attr = "signal") <- node.signal
      in.sign.impacts <- vector("list")
      signal.base.denoms <- vector("list")
      if (length(parent.nodes) > 1) {
        for (parent in parent.nodes) {
          if (!is.null(eval.exprs[[parent]])) {
            in.sign.impacts[[parent]] <- sprintf("(%s)^(1+I['%s','%s'])",
                                                paste("eval.exprs[['", parent, "']]",
                                                      sep = ""), parent, node)
            signal.base.denoms[[parent]] <- sprintf("(%s)",
                                                   paste("eval.exprs[['", parent, "']]",
                                                         sep = ""))
          }
          else {
            in.sign.impacts[[parent]] <- sprintf("(E['%s',1])^(1+I['%s','%s'])",
                                                parent, parent, node)
            signal.base.denoms[[parent]] <- sprintf("(E['%s',1])",
                                                   parent)
          }
        }
        signal.base.denom <- paste(signal.base.denoms,
                                  collapse = "+")
        signal.base <- sprintf("E[%d,1]/(%s)", i, signal.base.denom)
        in.signal.impact <- paste(in.sign.impacts, collapse = "+")
        eval.exprs[[node]] <- sprintf("(%s)*(%s)", signal.base,
                                     in.signal.impact)
      }
      else {
        parent <- parent.nodes[1]
        if (!is.null(eval.exprs[[parent]])) {
          in.signal.impact <- sprintf("(%s)^(I['%s','%s'])",
                                     paste("eval.exprs[['", parent, "']]", sep = ""),
                                     parent, node)
        }
        else {
          in.signal.impact <- sprintf("(E['%s',1])^(I['%s','%s'])",
                                     parent, parent, node)
        }
        signal.base <- sprintf("E[%d,1]", i)
        eval.exprs[[node]] <- sprintf("(%s)*(%s)", signal.base,
                                     in.signal.impact)
      }
    }
    else {
      eval.exprs[[node]] <- sprintf("E[%d,1]", i)
    }
  }

  signal.at.nodes <- unlist(nodeData(g, attr="signal"))
  names(signal.at.nodes) <- nodeIDs
  signal.at.sinks <- NULL
  if (!is.null(sink.nodes) && length(sink.nodes)>0 )
  {
    signal.at.sinks <- unlist(nodeData(g, sink.nodes, attr = "signal"))
    names(signal.at.sinks) <- nodeIDs[names(signal.at.sinks)]
  }

  return( list( signal.at.nodes=signal.at.nodes, signal.at.sinks=signal.at.sinks ) )
}







pipeline.PSFcalculation <- function(env)
{
  data(kegg.collection)
  
  if( is.null(env$indata.ensID.m) )
  {
    env$indata.ensID.m <- env$indata[env$gene.info$ensembl.mapping[,1],]
    env$indata.ensID.m <- do.call(rbind, by(env$indata.ensID.m, env$gene.info$ensembl.mapping[,2], colMeans))
  }

  if( ncol(env$indata) <= 100 )
  {
    util.info("Performing PSF calculation for samples:" )
    progressbar <- newProgressBar(min = 0, max = length(kegg.collection)); cat("\r")
    
    env$psf.results.samples <- list()

    for (i in 1:length(kegg.collection))
    {
      env$psf.results.samples[[names(kegg.collection)[i]]] <- list()

      for( m in 1:ncol(env$indata.ensID.m) )
      {
        if (!length(kegg.collection[[i]]) == 0)
        {
          g <- pathway.expression.mapping( kegg.collection[[i]]$graph, fc.m = 10^env$indata.ensID.m[,m] )

          env$psf.results.samples[[names(kegg.collection)[i]]][[colnames(env$indata.ensID.m)[m]]] <- psf.flow(g, node.ordering=kegg.collection[[i]]$order, sink.nodes=kegg.collection[[i]]$sink.nodes, nodeIDs=kegg.collection[[i]]$nodeIDs )
        }
      }

      setTxtProgressBar( progressbar, progressbar$getVal()+1 )
    }
    progressbar$kill()
  }



  if(length(unique(env$group.labels))>1)
  {
    util.info("Performing PSF calculation for groups:" )
    
    progressbar <- newProgressBar(min = 0, max = length(kegg.collection)); cat("\r")

    indata.ensID.groups <- do.call(cbind, by(t(env$indata.ensID.m),env$group.labels,colMeans)[unique(env$group.labels)])

    env$psf.results.groups <- list()

    for (i in 1:length(kegg.collection))
    {
      env$psf.results.groups[[names(kegg.collection)[i]]] <- list()

      for( m in 1:ncol(indata.ensID.groups) )
      {
        if (!length(kegg.collection[[i]]) == 0)
        {
          g <- pathway.expression.mapping( kegg.collection[[i]]$graph, fc.m = 10^indata.ensID.groups[,m] )

          env$psf.results.groups[[names(kegg.collection)[i]]][[colnames(indata.ensID.groups)[m]]] <- psf.flow(g, node.ordering=kegg.collection[[i]]$order, sink.nodes=as.character(kegg.collection[[i]]$sink.nodes), nodeIDs=kegg.collection[[i]]$nodeIDs )
        }
      }

      setTxtProgressBar( progressbar, progressbar$getVal()+1 )
    }
    progressbar$kill()
  }
  
  return(env)
}



pipeline.PSFoutput <- function(env)
{
  data(kegg.collection)
  
  if( !is.null(env$psf.results.samples) || !is.null(env$psf.results.groups) )
  {
    dir.create(file.path(paste(env$files.name, "- Results"),"PSF Analysis"), showWarnings=FALSE)
  }
  
  if( !is.null(env$psf.results.samples) )
  {
    output.path <- file.path(paste(env$files.name, "- Results"),"PSF Analysis","Sample Centered")
    
    dir.create(output.path, showWarnings=FALSE)
    
    psf.overview.heatmaps(env$psf.results.samples, output.path, env$group.colors, env$color.palette.heatmaps)
    psf.report.sheets(env,env$psf.results.samples, output.path, env$group.colors)
        
  }
  
  
  if( !is.null(env$psf.results.groups) )
  {
    output.path <- file.path(paste(env$files.name, "- Results"),"PSF Analysis","Group Centered")
    
    dir.create(output.path, showWarnings=FALSE)
    
    psf.overview.heatmaps(env$psf.results.groups, output.path, env$groupwise.group.colors, env$color.palette.heatmaps)
    psf.report.sheets(env,env$psf.results.groups, output.path, env$groupwise.group.colors)    
  }  
  
}
