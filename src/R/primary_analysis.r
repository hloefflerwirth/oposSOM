pipeline.prepareIndata <- function()
{
  indata.sample.mean <<- colMeans(indata)

  if (preferences$sample.quantile.normalization)
  {
    indata <<- Quantile.Normalization(indata)
  }

  colnames(indata) <<- make.unique(colnames(indata))
  names(group.labels) <<- make.unique(names(group.labels))
  names(group.colors) <<- make.unique(names(group.colors))


  indata.gene.mean <<- rowMeans(indata)

  if (preferences$feature.centralization)
  {
    indata <<- indata - indata.gene.mean
  }
}


pipeline.generateSOM <- function()
{
  som.result <<- som.linear.init(indata,somSize=preferences$dim.1stLvlSom)
  
  # Rotate/Flip First lvl SOMs

  if (preferences$rotate.SOM.portraits > 0)
  {
    for (i in 1:preferences$rotate.SOM.portraits)
    {
      o <- matrix(c(1:(preferences$dim.1stLvlSom^2)), preferences$dim.1stLvlSom, preferences$dim.1stLvlSom, byrow=TRUE)
      o <- o[rev(1:preferences$dim.1stLvlSom),]
      som.result <<- som.result[as.vector(o),]
    }
  }

  if (preferences$flip.SOM.portraits)
  {
    o <- matrix(c(1:(preferences$dim.1stLvlSom^2)), preferences$dim.1stLvlSom, preferences$dim.1stLvlSom, byrow=TRUE)
    som.result <<- som.result[as.vector(o),]
  }

  som.result <<- som.training( indata, som.result, prolongationFactor = preferences$training.extension, verbose = TRUE )
    
  metadata <<- som.result$weightMatrix
  colnames(metadata) <<- colnames(indata)

  som.result$weightMatrix <<- NULL


  ## set up SOM dependent variables
  
  gene.info$coordinates <<- apply( som.result$node.summary[som.result$feature.BMU,c("x","y")], 1, paste, collapse=" x " )
  names(gene.info$coordinates) <<- rownames(indata)
  
}
