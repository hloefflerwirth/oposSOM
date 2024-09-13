#' Calculates PSF formulas for each node in graphNEL object.
#' @param g graphNEL pathway object.
#' @param node.ordering order of nodes calculated with order.nodes function.
#' @param sink.nodes list of terminal (sink) nodes calculated with determine.sink.nodes function.
#' @param split logical, if true then the incoming signal will be proportionally splitted among the edges.
#' @param sum logical, default value is FALSE. When set to true pathway activity formulas will be calculated via addition, when set to false then activity formulas will be calculated via multiplication.
#' @param tmm_mode when set to true specific PSF configuration will be used for calculation of the pathway activity formulas described in https://www.frontiersin.org/articles/10.3389/fgene.2021.662464/full
#' @param tmm_update_mode when set to true specific PSF configuration will be used for calculation the pathway activity formulas described in https://www.frontiersin.org/articles/10.3389/fgene.2021.662464/full
#' @import graph

eval_formulas <- function(g, node.ordering, sink.nodes, split = TRUE, sum = FALSE, tmm_mode = FALSE, tmm_update_mode = FALSE) {
  
  node.order <- node.ordering$node.order
  node.rank <- node.ordering$node.rank
  
  nods <- names(node.order)
  eval.exprs = vector("list")
  
  #expressions
  node_fun = data.frame(as.character(graph::nodeData(g, nods, attr = "psf_function")), row.names=nods)
  
  #impacts
  I = matrix(data=NA, nrow=length(nods),ncol= length(nods))
  W = matrix(data = NA, nrow = length(nods), ncol = length(nods))
  rownames(I) = nods
  colnames(I) = nods
  rownames(W) = nods
  colnames(W) = nods
  
  for (node in nods){
    #     show(node)
    l = length(graph::edgeL(g)[[node]]$edges)
    for (e in 1:l){
      from = node
      to = graph::nodes(g)[graph::edgeL(g)[node][[1]]$edges][e]
      if (!is.na(to) && length(to)>0){
        weight = graph::edgeData(g, from = from, to =to, attr = "weight")[[1]]
        W[from, to] = weight
        impact = graph::edgeData(g, from = from, to =to, attr = "impact")[[1]]
        I[from, to] = impact
      }
    }
  }
  
  recalc.sinks =F
  if(is.null(sink.nodes)){
    recalc.sinks = T
    cat("\nsink nodes were not supplied. Those will be recalculated ...\n")
  }
  
  for (i in 1:length(nods)){
    
    parent.nodes <- graph::inEdges(nods[i], g)[[1]]
    
    # if(tmm_mode) {
    #   ## skipping nodes with FC value 1
    #   parent.nodes <- parent.nodes[which(unlist(graph::nodeData(g, parent.nodes, attr = "signal")) != 1)]
    # }
    
    node = nods[i]
    
    if (length(parent.nodes)>0) {
      #       node.exp <- nodeData(g, i, attr = "expression")[[1]]
      # node.exp = E[i,1]
      
      ### input signal processing of function nodes (will be implemented in the future)
      # if("psf_function" %in% names(graph::nodeDataDefaults(g))) {
      #   if(unname(unlist(graph::nodeData(g, node, attr = "psf_function"))) %in% c("min", "max", "sum")) {
      #     in.signal <- get(unname(unlist(graph::nodeData(g, node, attr = "psf_function"))), envir = globalenv())(unlist(graph::nodeData(g, parent.nodes, attr = "signal")))
      #   } else {
      #     in.signal <- unlist(graph::nodeData(g, parent.nodes, attr = "signal"))
      #   }
      # } else {
      #   in.signal <- unlist(graph::nodeData(g, parent.nodes, attr = "signal"))
      # }
      
      if(sum){
        
        in.sign.impacts = vector("list")
        for (parent in parent.nodes){
          if(!is.null(eval.exprs[[parent]])){
            in.sign.impacts[[parent]] = sprintf("(W['%s','%s'])*(%s)*(I['%s','%s'])", parent, node, paste("eval.exprs[['",parent,"']]",sep=""), parent, node)
          } else {
            in.sign.impacts[[parent]] = sprintf("(W['%s','%s'])*(E['%s',1])*(I['%s','%s'])", parent, node, parent, parent, node)
          }
        }
        in.signal.impact = paste(in.sign.impacts,collapse="+")
        
        eval.exprs[[node]] = sprintf("(E[%d,1])+(%s)", i, in.signal.impact)
        
        
      } else {
        if(split){
          
          in.sign.impacts = vector("list")
          signal.base.denoms = vector("list")
          if(length(parent.nodes) > 1) {
            for (parent in parent.nodes){
              if(!is.null(eval.exprs[[parent]])){
                in.sign.impacts[[parent]] = sprintf("(W['%s','%s'])*(%s)^(1+I['%s','%s'])", parent, node, paste("eval.exprs[['",parent,"']]",sep=""), parent, node)
                signal.base.denoms[[parent]] = sprintf("(%s)",paste("eval.exprs[['",parent,"']]",sep=""))
              } else {
                in.sign.impacts[[parent]] = sprintf("(W['%s','%s'])*(E['%s',1])^(1+I['%s','%s'])", parent, node, parent, parent, node)
                signal.base.denoms[[parent]] = sprintf("(E['%s',1])",parent)
              }
            }
            signal.base.denom = paste(signal.base.denoms, collapse="+")
            signal.base = sprintf("E[%d,1]/(%s)", i, signal.base.denom)
            in.signal.impact = paste(in.sign.impacts,collapse="+")
            eval.exprs[[node]] = sprintf("(%s)*(%s)", signal.base, in.signal.impact)
          } else {
            parent = parent.nodes[1]
            if(!is.null(eval.exprs[[parent]])){
              in.signal.impact = sprintf("(W['%s','%s'])*(%s)^(I['%s','%s'])", parent, node, paste("eval.exprs[['",parent,"']]",sep=""), parent, node)
            } else {
              in.signal.impact = sprintf("(W['%s','%s'])*(E['%s',1])^(I['%s','%s'])", parent, node, parent, parent, node)
            }
            signal.base = sprintf("E[%d,1]", i)
            eval.exprs[[node]] = sprintf("(%s)*(%s)", signal.base, in.signal.impact)
          }
          
          #### special formula for signal exponential decay, not used yet ####
          # a = 2000
          # proportion = in.signal/sum(in.signal)
          # node.signal <- sum((proportion * weight * node.exp) *  a*(2/(1+exp((-2*in.signal^impact)/a)) - 1))
          
          
        } else {
          
          in.sign.impacts = vector("list")
          signal.base.denoms = vector("list")
          
          if(tmm_mode) {
            
            for (parent in parent.nodes){
              if(!is.null(eval.exprs[[parent]])){
                in.sign.impacts[[parent]] = sprintf("(W['%s','%s'])*(%s)^(I['%s','%s'])", parent, node, paste("eval.exprs[['",parent,"']]",sep=""), parent, node)
                signal.base.denoms[[parent]] = sprintf("(%s)",paste("eval.exprs[['",parent,"']]",sep=""))
              } else {
                in.sign.impacts[[parent]] = sprintf("(W['%s','%s'])*(E['%s',1])^(I['%s','%s'])", parent, node, parent, parent, node)
                signal.base.denoms[[parent]] = sprintf("(E['%s',1])",parent)
              }
            }
            
            if(node_fun[i,1] %in% c("min", "max", "sum")) {
              signal.base.denom <- paste(signal.base.denoms, collapse = ", ")
              
              in.signal <- sprintf("%s(%s, %s)", node_fun[i,1], signal.base.denom, "na.rm = TRUE")
              
              eval.exprs[[node]] <- sprintf("(E[%d,1])*(%s)", i, in.signal)
              
            } else {
              
              if(tmm_update_mode) {
                total_signal <- paste(c(sprintf("E[%d,1]", i), in.sign.impacts), collapse = ",")
                
                eval.exprs[[node]] <- sprintf("%s(%s, %s)", "prod", total_signal, "na.rm = TRUE")
                
                # eval.exprs[[node]] <- paste(c(sprintf("E[%d,1]", i), in.sign.impacts), collapse = "*")
              } else {
                total_signal <- paste("(", paste(sprintf("E[%d,1]", i), in.sign.impacts, sep = "*"), ")", collapse = ", ")
                
                eval.exprs[[node]] <- sprintf("%s(%s, %s)", "prod", total_signal, "na.rm = TRUE")
                # eval.exprs[[node]] <- paste("(", paste(sprintf("E[%d,1]", i), in.sign.impacts, sep = "*"), ")", collapse = "*")
                
              }
              
            }
            
          } else {
            #Returns the product of signals without splitting - is for updateing, but applies only in this special case where all the rules are s*t, or 1/s*t
            
            in.sign.impacts = vector("list")
            if(length(parent.nodes) > 1) {
              for (parent in parent.nodes){
                if(!is.null(eval.exprs[[parent]])){
                  in.sign.impacts[[parent]] = sprintf("(W['%s','%s'])*(%s)^(I['%s','%s'])", parent, node, paste("eval.exprs[['",parent,"']]",sep=""), parent, node)
                } else {
                  in.sign.impacts[[parent]] = sprintf("(W['%s','%s'])*(E['%s',1])^(I['%s','%s'])", parent, node, parent, parent, node)
                }
              }
              
              in.signal.impact = paste(in.sign.impacts,collapse="*")
              signal.base = sprintf("E[%d,1]", i)
              
              eval.exprs[[node]] = sprintf("(%s)*(%s)", signal.base, in.signal.impact)
            } else {
              parent = parent.nodes[1]
              if(!is.null(eval.exprs[[parent]])){
                in.signal.impact = sprintf("(W['%s','%s'])*(%s)^(I['%s','%s'])", parent, node, paste("eval.exprs[['",parent,"']]",sep=""), parent, node)
              } else {
                in.signal.impact = sprintf("(W['%s','%s'])*(E['%s',1])^(I['%s','%s'])", parent, node, parent, parent, node)
              }
              signal.base = sprintf("E[%d,1]", i)
              eval.exprs[[node]] = sprintf("(%s)*(%s)", signal.base, in.signal.impact)
            }
            
          }
          
          
        }
      }
      
    } else {
      eval.exprs[[node]] = sprintf("E[%d,1]",i)
    }
    
    if(recalc.sinks){
      if (length(parent.nodes) > 0){
        child.nodes <- unlist(graph::edges(g, names(node.order)[i]))
        if (length(child.nodes) >0) {
          child.node.ranks <- node.rank[child.nodes]
          rank.comp = child.node.ranks > node.rank[i]
          if (all(rank.comp))
            sink.nodes <- c(sink.nodes, names(node.order)[i])
        } else {
          sink.nodes <- c(sink.nodes, names(node.order)[i])
        }
      }
    }
    # sink.nodes <- c(sink.nodes, names(node.order)[i])
    
  }
  
  return(list("graph" = g, "order"=node.ordering, "sink.nodes" = sink.nodes,
              "eval.exprs" = eval.exprs, "I" = I, "W" = W))
  
}

### running psf on pathway

psf_processing <- function(fc.expression, pathway) 
{
  
  ### building psf formulas
  eval_g <- eval_formulas(g = pathway$graph, node.ordering = pathway$order, sink.nodes = pathway$sink.nodes, split = T, sum = F, tmm_mode = F, tmm_update_mode = F)
  
  
  ### building exp FC matrix for the pathway
  gene.data <- graph::nodeData(pathway$graph)[names(pathway$order$node.order)]
  pathway_exp_data <- lapply(gene.data, function(x) {
    
    genes.in.node = which(rownames(fc.expression) %in% x$genes)
    expression.values <- fc.expression[genes.in.node,, drop = F]
    
    if (nrow(expression.values)>0){
      if (x$type == "gene"){
        expression <- matrix(colMeans(expression.values), ncol = ncol(fc.expression))
      } else { 
        if(x$type == "group") {
          expression <- apply(expression.values, 2, FUN = min)
        } else {
          expression <- matrix(rep(1, ncol(fc.expression)), ncol = ncol(fc.expression))
        }
      }
    } else {
      expressions <- matrix(rep(1, ncol(fc.expression)), ncol = ncol(fc.expression))
    }
    
  })
  
  pathway_exp_data <- Reduce(rbind, pathway_exp_data)
  
  rownames(pathway_exp_data) <- names(pathway$order$node.order)
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
    
  
  # pathway$psf_activities <- psf_activities
  # pathway$exp_fc <- pathway_exp_data
  # 
  # 
  # pathway$eval.exprs <- eval_g$eval.exprs
  # 
  # pathway$I <- I
  # pathway$W <- W
  
  # pathway$graph <- NULL
  
  return(psf_activities)
  
}



##### pipeline functions ####

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
