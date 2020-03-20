newProgressBar <- function (min = 0, max = 1, initial = 0) 
{
  .val <- initial
  .nb <- 0L
  .pc <- -1L
  width <- 70

  up3 <- function(value) {
    if (!is.finite(value) || value < min || value > max) 
      return()
    .val <<- value
    nb <- round(width * (value - min)/(max - min))
    pc <- round(100 * (value - min)/(max - min))
    if (nb == .nb && pc == .pc) 
      return()
    cat(paste0("\r  |", rep(" ", 1 * width + 6)))
    cat(paste(c("\r  |", rep.int("=", nb), rep.int(" ", 1 * (width - nb)), sprintf("| %3d%%", pc)), collapse = "") )
    flush.console()
    .nb <<- nb
    .pc <<- pc
  }
  getVal <- function() .val
  kill <- function() { cat("\n");  flush.console() }
  structure(list(getVal = getVal, up = up3, kill=kill), class = "txtProgressBar")
}



som.training.phase <- function(indata, weightMatrix, metricSamples, epochs, 
                               initLearnRate, inverseLearnRate, initRadius, progressbar )
{
  cycleIntern <- 0
  maxCycleIntern <- epochs * nrow(indata)
  
  deltaMatrix <- matrix(0, nrow=nrow(weightMatrix), ncol=ncol(weightMatrix))
  distanceVector  <- rep(0, nrow(weightMatrix))
  neighborhoodVector  <- rep(0, nrow(weightMatrix))

  somSize <- sqrt(nrow(weightMatrix))
  
  any.na <- any(is.na(indata))

  for (epoch.i in seq(epochs-1) )
  { 
    t1 <- Sys.time()
    for (it in seq(nrow(indata)) )
    {
      cycleIntern <- cycleIntern + 1
      radius.t <- 1+(initRadius-1)*(maxCycleIntern-cycleIntern)/(maxCycleIntern)
      learnRate.t <- initLearnRate * inverseLearnRate / (inverseLearnRate + cycleIntern)

      calculateDelta( weightMatrix, indata[it,], any.na, deltaMatrix )
      calculateEuclideanDistances( deltaMatrix[,metricSamples], distanceVector )
       
      BMU <- which.min(distanceVector)
      calculateNeighborhoodMatrix( BMU, somSize, radius.t, neighborhoodVector )
      
      weightMatrix <- weightMatrix + learnRate.t * (deltaMatrix * neighborhoodVector)
    }
      
    if(!is.null(progressbar))
    {
      if( progressbar$getVal() == 0 )
      {
        past.runtime = as.double(difftime( Sys.time(),t1, units="secs"))
        util.info("Remaining time for SOM training: ~", ceiling(11*past.runtime/60), "min = ~", round(11*past.runtime/3600,1),"h")      
      } 
      setTxtProgressBar( progressbar, progressbar$getVal()+1 )
    }
  }
  

  # remember BMU in last epoch
  BMU <- rep(NA,nrow(indata))
  names(BMU) <- rownames(indata)
  for (it in seq(nrow(indata)) )
  {
    cycleIntern <- cycleIntern + 1
    radius.t <- 1+(initRadius-1)*(maxCycleIntern-cycleIntern)/(maxCycleIntern)
    learnRate.t <- initLearnRate * inverseLearnRate / (inverseLearnRate + cycleIntern)
    
    calculateDelta( weightMatrix, indata[it,], any.na, deltaMatrix )
    calculateEuclideanDistances( deltaMatrix[,metricSamples], distanceVector )
    
    BMU[it] <- which.min(distanceVector)
    calculateNeighborhoodMatrix( BMU[it], somSize, radius.t, neighborhoodVector )
    
    weightMatrix <- weightMatrix + learnRate.t * (deltaMatrix * neighborhoodVector)
  }
  if(!is.null(progressbar)) setTxtProgressBar( progressbar, progressbar$getVal()+1 )
  
  return( list( weightMatrix=weightMatrix, BMU=BMU ) )
}





som.training <- function( indata, weightMatrix, metricSamples, prolongationFactor = 1, verbose = FALSE )
{
  if( missing(metricSamples) ) metricSamples <- seq( ncol(indata) )

  if(verbose) { pb <-newProgressBar(min = 0, max = 12); cat("\r") } else pb <- NULL
  
  somSize <- sqrt(nrow(weightMatrix))

  som.result <- som.training.phase( indata, weightMatrix, metricSamples,
                              epochs=2*prolongationFactor,
                              initLearnRate=0.05,
                              inverseLearnRate=nrow(indata)*2*prolongationFactor/100,
                              initRadius=somSize, pb )

  som.result <- som.training.phase( indata, som.result$weightMatrix, metricSamples,
                              epochs=10*prolongationFactor,
                              initLearnRate=0.02,
                              inverseLearnRate=nrow(indata)*10*prolongationFactor/100,
                              initRadius=min(3,somSize), pb )

  node.summary <- data.frame(x = rep(1:somSize, somSize), 
                             y = rep(1:somSize, each=somSize),
                             n.features = NA )
  node.summary$n.features <- as.vector( table( som.result$BMU )[as.character(1:(somSize^2))] )
  node.summary$n.features[ which(is.na(node.summary$n.features)) ] <- 0
  
  if(verbose) pb$kill()

  return(  list(
      "weightMatrix" = som.result$weightMatrix,
      "node.summary" = node.summary,
      "feature.BMU" = som.result$BMU
    ) )
}

