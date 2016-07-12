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
  som.result <<- som.init(indata, xdim=preferences$dim.1stLvlSom, ydim=preferences$dim.1stLvlSom, init="linear")

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


  # Train SOM

  # The following would train the SOM in one step...
  #som.result <<- som(indata, xdim=preferences$dim.1stLvlSom, ydim=preferences$dim.1stLvlSom)

  # We split training in two steps to estimate the time we need.
  t1 <- system.time({
    som.result <<- som.train(indata, som.result, xdim=preferences$dim.1stLvlSom,
                             ydim=preferences$dim.1stLvlSom, alpha=0.05,
                             radius=preferences$dim.1stLvlSom,
                             rlen=nrow(indata)*2*preferences$training.extension,
                             inv.alp.c=nrow(indata)*2*preferences$training.extension/100)
  })

  util.info("Remaining time for SOM training: ~", ceiling(5*t1[3]/60), "min = ~", round(5*t1[3]/3600,1),"h")

  som.result <<- som.train(indata, som.result$code, xdim=preferences$dim.1stLvlSom,
                           ydim=preferences$dim.1stLvlSom, alpha=0.02,
                           radius=min(3, preferences$dim.1stLvlSom),
                           rlen=nrow(indata)*10*preferences$training.extension,
                           inv.alp.c=nrow(indata)*10*preferences$training.extension/100)

  metadata <<- som.result$code
  colnames(metadata) <<- colnames(indata)

  som.result$data <<- NULL
  som.result$code <<- NULL


  ## set up SOM dependent variables
  gene.info$coordinates <<- apply(som.result$visual[,c(1,2)]+1, 1, paste, collapse=" x ")
  names(gene.info$coordinates) <<- rownames(indata)

  som.result$nodes <<- (som.result$visual[,"x"] + 1) + som.result$visual[,"y"] * preferences$dim.1stLvlSom
  names(som.result$nodes) <<- rownames(indata)
}
