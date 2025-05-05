
eval_formulas <- function( pathway ) {
  
  ordered.nodes <- names(pathway$node.order)
  eval.exprs <- vector("list")
  
  #impacts
  I <- matrix( NA, nrow=length(ordered.nodes),ncol=length(ordered.nodes), dimnames=list(ordered.nodes,ordered.nodes) )
  W <- matrix( NA, nrow=length(ordered.nodes),ncol=length(ordered.nodes), dimnames=list(ordered.nodes,ordered.nodes) )

  for( edge in pathway$edge.info )  
  {
    W[ edge$entry1, edge$entry2 ] <- edge$weight
    I[ edge$entry1, edge$entry2 ] <- edge$impact
  }
  
  for( node.i in seq(ordered.nodes) )
  {
    node.id <- ordered.nodes[ node.i ]
    parent.nodes <- sapply(pathway$edge.info,"[[","entry1")[ which( sapply(pathway$edge.info,"[[","entry2")==node.id ) ]
      
    if( length(parent.nodes)>0 ) 
    {
      # incoming signal will be proportionally splitted among the edges
      
      in.sign.impacts <- vector("list")
      signal.base.denoms <- vector("list")
      
      if(length(parent.nodes) > 1) 
      {
        for (parent in parent.nodes)
        {
          if(!is.null(eval.exprs[[parent]]))
          {
            in.sign.impacts[[parent]] = sprintf("(W['%s','%s'])*(%s)^(1+I['%s','%s'])", parent, node.id, paste("eval.exprs[['",parent,"']]",sep=""), parent, node.id)
            signal.base.denoms[[parent]] = sprintf("(%s)",paste("eval.exprs[['",parent,"']]",sep=""))
            
          } else 
          {
            in.sign.impacts[[parent]] = sprintf("(W['%s','%s'])*(E['%s',1])^(1+I['%s','%s'])", parent, node.id, parent, parent, node.id)
            signal.base.denoms[[parent]] = sprintf("(E['%s',1])",parent)
          }
        }
        signal.base.denom <- paste(signal.base.denoms, collapse="+")
        signal.base <- sprintf("E[%d,1]/(%s)", node.i, signal.base.denom)
        in.signal.impact = paste(in.sign.impacts,collapse="+")
        eval.exprs[[node.id]] = sprintf("(%s)*(%s)", signal.base, in.signal.impact)
        
      } else 
      {
        if(!is.null(eval.exprs[[parent.nodes]])){
          in.signal.impact <- sprintf("(W['%s','%s'])*(%s)^(I['%s','%s'])", parent.nodes, node.id, paste("eval.exprs[['",parent.nodes,"']]",sep=""), parent.nodes, node.id)
        } else {
          in.signal.impact <- sprintf("(W['%s','%s'])*(E['%s',1])^(I['%s','%s'])", parent.nodes, node.id, parent.nodes, parent.nodes, node.id)
        }
        signal.base <- sprintf("E[%d,1]", node.i)
        eval.exprs[[node.id]] = sprintf("(%s)*(%s)", signal.base, in.signal.impact)
      }
      
    } else 
    {
      eval.exprs[[node.id]] <- sprintf("E[%d,1]",node.i)
    }
        
  }
  
  return(list("eval.exprs"=eval.exprs, "I"=I, "W"=W))
  
}

### running psf on pathway

psf_processing <- function(fc.expression, pathway) 
{
  
  ### building psf formulas
  eval_g <- eval_formulas( pathway )
  
  ### building exp FC matrix for the pathway
  get.gene.node.expression <- function(ens.ids,fc.expression)
  {
    genes.in.node <- intersect(rownames(fc.expression), ens.ids)
    
    if( length(genes.in.node)>0 )
    {
      expression.values <- fc.expression[genes.in.node,, drop = F]
      return( matrix( apply(expression.values, 2, max), ncol = ncol(fc.expression)) )
      
    } else    
    return( matrix(rep(1, ncol(fc.expression)), ncol = ncol(fc.expression)) )
  }
  
  pathway_exp_data <- lapply( pathway$node.info[names(pathway$node.order)], function(x) 
  {
    expression <- matrix(rep(1, ncol(fc.expression)), ncol = ncol(fc.expression))
    
    if (x$type == "gene")
    {
      expression <- get.gene.node.expression( x$ens.ids, fc.expression )

    } else 
    if(x$type == "group") 
    {
      gr.expression <- lapply( x$components, function(comp) get.gene.node.expression( comp$ens.ids, fc.expression ) )
      gr.expression <- do.call( rbind, gr.expression )  
      gr.expression <- gr.expression[ which(apply(gr.expression,1,function(y)any(y!=1))), , drop=F]
      
      if(nrow(gr.expression)>0)
        expression <- apply(gr.expression, 2, min)
    } 
    
    return(expression)    
  })
  
  pathway_exp_data <- Reduce(rbind, pathway_exp_data)
  
  rownames(pathway_exp_data) <- names(pathway$node.order)
  colnames(pathway_exp_data) <- colnames(fc.expression)
  

  ### calculating psf with pathway formulas
  I = eval_g$I
  W = eval_g$W
  eval.exprs <- eval_g$eval.exprs
  
  eval_calculator <- function(exp_ind, exp_mat) {
    E <- exp_mat[,exp_ind, drop = F]
    for (i in 1:length(eval.exprs)){
      eval.exprs[[i]] = eval(parse(text = eval.exprs[[i]]))
    }
    return(eval.exprs)
  }
  
  psf_activities <- sapply(1:ncol(pathway_exp_data), function(x) eval_calculator(exp_ind = x, exp_mat = pathway_exp_data) )
  
  psf_activities <- apply(psf_activities, 2, as.numeric)
  colnames(psf_activities) <- colnames(pathway_exp_data)
  rownames(psf_activities) <- rownames(pathway_exp_data)
    
  
  return(psf_activities)
}



##### pipeline functions ####

pipeline.PSFcalculation <- function(env)
{
  data(kegg.collection)
  
  all.psf.genes <- unique( do.call( c, lapply(kegg.collection, function(x) unique( do.call( c, sapply(x$node.info,function(xx)xx$ens.ids) ) ) ) ) )
  if( mean( all.psf.genes %in% env$gene.info$ids ) < 0.50 )
  {
    util.warn("Too few KEGG signalling pathway genes covered. Disabling PSF analysis.")
    env$preferences$activated.modules$psf.analysis <- FALSE
    return(env)
  }
  
  if( is.null(env$indata.ensID.m) )
  {
    env$indata.ensID.m <- env$indata[env$gene.info$ensembl.mapping[,1],]
    env$indata.ensID.m <- do.call(rbind, by.minicluster(env$indata.ensID.m, env$gene.info$ensembl.mapping[,2], colMeans))
  }
  
  if( ncol(env$indata) <= 100 )
  {
    util.info("Performing PSF calculation for samples:" )
    progressbar <- newProgressBar(min = 0, max = length(kegg.collection)); cat("\r")
    
    fc.expression <- 10^env$indata.ensID.m
    
    env$psf.results.samples <- list()
    
    for (i in 1:length(kegg.collection))
    {
      env$psf.results.samples[[names(kegg.collection)[i]]] <- psf_processing(fc.expression, kegg.collection[[i]]) 
      
      setTxtProgressBar( progressbar, progressbar$getVal()+1 )
    }
    progressbar$kill()
  }
  
  
  if(length(unique(env$group.labels))>1)
  {
    util.info("Performing PSF calculation for groups:" )
    
    progressbar <- newProgressBar(min = 0, max = length(kegg.collection)); cat("\r")
    
    indata.ensID.groups <- do.call(cbind, by(t(env$indata.ensID.m),env$group.labels,colMeans)[unique(env$group.labels)])
    
    fc.expression <- 10^indata.ensID.groups
    
    env$psf.results.groups <- list()
    
    for (i in 1:length(kegg.collection))
    {
      env$psf.results.groups[[names(kegg.collection)[i]]] <- psf_processing(fc.expression, pathway=kegg.collection[[i]]) 
      
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
    dir.create("PSF Analysis", showWarnings=FALSE)
  }
  
  if( !is.null(env$psf.results.samples) )
  {
    output.path <- file.path("PSF Analysis","Sample Centered")
    
    dir.create(output.path, showWarnings=FALSE)
    
    psf.overview.heatmaps(psf.results=env$psf.results.samples, output.path, group.colors=env$group.colors, color.palette=env$color.palette.heatmaps)
    psf.report.sheets(env,psf.results=env$psf.results.samples, output.path, bar.colors=env$group.colors)
  }
  
  
  if( !is.null(env$psf.results.groups) )
  {
    output.path <- file.path("PSF Analysis","Group Centered")
    
    dir.create(output.path, showWarnings=FALSE)
    
    psf.overview.heatmaps(env$psf.results.groups, output.path, group.colors=env$groupwise.group.colors, env$color.palette.heatmaps)
    psf.report.sheets(env,psf.results=env$psf.results.groups, output.path, bar.colors=env$groupwise.group.colors)
  }
  
}
