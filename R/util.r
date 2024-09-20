# outputs a message
util.cat <- function(..., suffix="\n") {
  message(paste(..., collapse=" "), suffix, appendLF=FALSE)
  flush.console()
}

# logs a message
util.log <- function(..., prefix, suffix="\n") {
  util.cat(c(paste("[", format(Sys.time()), "][", prefix, "]", sep=""), ...),
           suffix=suffix)
}

# logs an info
util.info <- function(..., suffix="\n") {
  util.log(c(...), prefix="INFO", suffix=suffix)
}

# logs a warning
util.warn <- function(..., suffix="\n") {
  util.log(c(...), prefix="WARN", suffix=suffix)
}

# logs an fatal error
util.fatal <- function(..., suffix="\n") {
  util.log(c(...), prefix="FATAL", suffix=suffix)
}

# outputs a progress bar
util.progress <- function(cur, max, len=48) {
  x <- 1

  if (max <= len) {
    x <- round((cur * len) / max)
  } else {
    x <- which(floor(seq(0, max, length.out=len+1)) == floor(cur))[1] - 1
  }

  if (is.na(x)) {
    return()
  }

  util.info(paste("[",
                  paste(rep("#", x), collapse=""),
                  paste(rep(" ", len - x), collapse=""),
                  "]", sep=""), suffix="\r")
}

# terminate progress bar
util.progress.terminate <- function(len=48) {
  util.info(paste("[", paste(rep("#", len), collapse=""), "]", sep=""))
}

#
util.call <- function(fn, env) {
  environment(fn) <- env
  return(fn())
}


by.minicluster <- function(data,ids,fun)
{
  ret <- lapply(unique(ids), function(x) fun( data[which(ids== x), , drop = FALSE] ) )
  names(ret) <- unique(ids)
  return(ret)
}          



chunk.apply.rows <- function(data,fun.para,add.data=NULL)
{
  if(env$preferences$max.cores<2) # run sequential apply
  {
    for(add.n in names(add.data) )
    {
      eval( parse(text = paste(add.n,"<- add.data[[add.n]]" ) ) )
    }
    
    environment(fun.para)<-environment()
    return( apply(data,1,fun.para) )
  }
  
  runcomp <- function()
  {
    plan(multisession, workers=min(env$preferences$max.cores,length(chunk.list$data)))
    on.exit(closeAllConnections())
    oopts <- options(future.globals.maxSize = 4.0 * 1e9)
    on.exit(options(oopts))
    
    chunk.results <- future_lapply( chunk.list$data, apply, 1, chunk.list$fun )
    
    if( is.matrix(chunk.results[[1]]) )
    {
      chunk.results <- do.call(rbind,chunk.results)
      colnames(chunk.results) <- colnames(data)
      
    } else
    {
      n <- do.call( c, sapply(chunk.results,names) )
      chunk.results <- do.call(c,chunk.results)
      names(chunk.results) <- n
    }
    return(chunk.results)
  }  
  
  # now following some crazy tricks with globals and variable evaluation to workaround slowdown of parallel computing
  
  assign("chunk.list",NULL,envir=.GlobalEnv)          
  fun.para <- paste( deparse(fun.para), collapse="\n")
  eval( parse(text = paste("fun.para <- ",fun.para)) )
  
  n.chunks <- env$preferences$max.cores
  chunk.list <<- split( rownames(data), cut(seq(nrow(data)), n.chunks, labels=F ) )
  chunk.list <<- lapply( chunk.list, function(x) data[x,] )
  chunk.list <<- list(fun=fun.para,data=chunk.list)
  
  for(add.n in names(add.data) )
  {
    eval( parse(text = paste(add.n,"<- add.data[[add.n]]" ) ) )
  }
  
  ret <- runcomp()
  
  rm(chunk.list,envir = .GlobalEnv)
  
  return(ret)
}





# chunk.apply.cols <- function(data,fun.para,add.data=NULL)
# {
#   if(env$preferences$max.cores<2) # run sequential apply
#   {
#     for(add.n in names(add.data) )
#     {
#       eval( parse(text = paste(add.n,"<- add.data[[add.n]]" ) ) )
#     }
#     
#     environment(fun.para)<-environment()
#     return( apply(data,2,fun.para) )
#   }
#   
#   runcomp <- function()
#   {
#     plan(multisession, workers=min(env$preferences$max.cores,length(chunk.list$data)))
#     on.exit(closeAllConnections())
#     oopts <- options(future.globals.maxSize = 4.0 * 1e9)
#     on.exit(options(oopts))
#     
#     chunk.results <- future_lapply( chunk.list$data, apply, 2, chunk.list$fun )
#     
#     if( is.matrix(chunk.results[[1]]) )
#     {
#       chunk.results <- do.call(cbind,chunk.results)
#       rownames(chunk.results) <- rownames(data)
#       
#     } else
#     {
#       n <- do.call( c, sapply(chunk.results,names) )
#       chunk.results <- do.call(c,chunk.results)
#       names(chunk.results) <- n
#     }
#     return(chunk.results)
#   }  
#   
#   # now following some crazy tricks with globals and variable evaluation to workaround slowdown of parallel computing
#   
#   assign("chunk.list",NULL,envir=.GlobalEnv)          
#   fun.para <- paste( deparse(fun.para), collapse="\n")
#   eval( parse(text = paste("fun.para <- ",fun.para)) )
#   
#   n.chunks <- env$preferences$max.cores
#   chunk.list <<- split( rownames(data), cut(seq(nrow(data)), n.chunks, labels=F ) )
#   chunk.list <<- lapply( chunk.list, function(x) data[x,] )
#   chunk.list <<- list(fun=fun.para,data=chunk.list)
#   
#   for(add.n in names(add.data) )
#   {
#     eval( parse(text = paste(add.n,"<- add.data[[add.n]]" ) ) )
#   }
#   
#   ret <- runcomp()
#   
#   rm(chunk.list,envir = .GlobalEnv)
#   
#   return(ret)
# }






chunk.sapply <- function(data,fun.para,add.data,max.chunks=env$preferences$max.cores)
{
  if(max.chunks<2 || env$preferences$max.cores<2) # run sequential sapply
  {
    for(add.n in names(add.data) )
    {
      eval( parse(text = paste(add.n,"<- add.data[[add.n]]" ) ) )
    }
    
    environment(fun.para)<-environment()
    return( sapply(data,fun.para) )
  }
  
  
  runcomp <- function()
  {
    plan(multisession, workers=min(env$preferences$max.cores,length(chunk.list$data)))
    on.exit(closeAllConnections())
    oopts <- options(future.globals.maxSize = 4.0 * 1e9)
    on.exit(options(oopts))
    
    chunk.results <- future_lapply( chunk.list$data, sapply, chunk.list$fun )
    
    if( is.matrix(chunk.results[[1]]) )
    {
      chunk.results <- do.call(cbind,chunk.results)
      
    } else
    {
      n <- do.call( c, sapply(chunk.results,names) )
      chunk.results <- do.call(c,chunk.results)
      names(chunk.results) <- n
    }
    
    return(chunk.results)
  }  
  
  # now following some crazy tricks with globals and variable evaluation to workaround slowdown of parallel computing
  
  assign("chunk.list",NULL,envir=.GlobalEnv)          
  fun.para <- paste( deparse(fun.para), collapse="\n")
  eval( parse(text = paste("fun.para <- ",fun.para)) )
  
  n.chunks <- min( env$preferences$max.cores, max.chunks )
  chunk.list <<- split( names(data), cut(seq(length(data)), n.chunks, labels=F ) )
  chunk.list <<- lapply( chunk.list, function(x) data[x] )
  chunk.list <<- list(fun=fun.para,data=chunk.list)
  
  for(add.n in names(add.data) )
  {
    eval( parse(text = paste(add.n,"<- add.data[[add.n]]" ) ) )
  }
  
  ret <- runcomp()
  
  rm(chunk.list,envir = .GlobalEnv)
  
  return(ret)
}

















