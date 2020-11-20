pipeline.prepareIndata <- function(env)
{
  env$indata.sample.mean <- colMeans(env$indata)

  if (env$preferences$sample.quantile.normalization)
  {
    env$indata <- Quantile.Normalization(env$indata)
  }

  colnames(env$indata) <- make.unique(colnames(env$indata))
  names(env$group.labels) <- make.unique(names(env$group.labels))
  names(env$group.colors) <- make.unique(names(env$group.colors))


  env$indata.gene.mean <- rowMeans(env$indata)

  if (env$preferences$feature.centralization)
  {
    env$indata <- env$indata - env$indata.gene.mean
  }
  return(env)
}


pipeline.generateSOM <- function(env)
{
  env$som.result <- som.linear.init(env$indata,somSize=env$preferences$dim.1stLvlSom)
  
  # Rotate/Flip First lvl SOMs

  if (env$preferences$rotate.SOM.portraits > 0)
  {
    for (i in 1:env$preferences$rotate.SOM.portraits)
    {
      o <- matrix(c(1:(env$preferences$dim.1stLvlSom^2)), env$preferences$dim.1stLvlSom, env$preferences$dim.1stLvlSom, byrow=TRUE)
      o <- o[rev(1:env$preferences$dim.1stLvlSom),]
      env$som.result <- env$som.result[as.vector(o),]
    }
  }

  if (env$preferences$flip.SOM.portraits)
  {
    o <- matrix(c(1:(env$preferences$dim.1stLvlSom^2)), env$preferences$dim.1stLvlSom, env$preferences$dim.1stLvlSom, byrow=TRUE)
    env$som.result <- env$som.result[as.vector(o),]
  }

  env$som.result <- som.training( env$indata, env$som.result, prolongationFactor = env$preferences$training.extension, verbose = TRUE )
    
  env$metadata <- env$som.result$weightMatrix
  colnames(env$metadata) <- colnames(env$indata)

  env$som.result$weightMatrix <- NULL


  ## set up SOM dependent variables
  
  env$gene.info$coordinates <- apply( env$som.result$node.summary[env$som.result$feature.BMU,c("x","y")], 1, paste, collapse=" x " )
  names(env$gene.info$coordinates) <- rownames(env$indata)
  
  return(env)
}
