som.linear.init.subdata <- function (indata, somSize) 
{
  if( !env$preferences$activated.modules$largedata.mode )
  {
    pca <- prcomp(indata)
    loadings.x <- pca$sdev[1] * pca$rotation[, 1]
    loadings.y <- pca$sdev[2] * pca$rotation[, 2]
    
  } else
  {
    indata.chunks <- split( colnames(indata), cut(seq(ncol(indata)), ncol(indata)/10, labels=F ) )
    condensed.indata <- do.call(cbind, lapply( indata.chunks, function(x) rowMeans(indata[,x]) )   )
      
    pca <- prcomp(condensed.indata)
    loadings.condensed.x <- pca$sdev[1] * pca$rotation[, 1]
    loadings.condensed.y <- pca$sdev[2] * pca$rotation[, 2]
    
    loadings.x <- do.call( c, lapply( names(loadings.condensed.x), function(chunk)
    {
      ret <- rep(loadings.condensed.x[chunk], length(indata.chunks[[chunk]]) )
      names(ret) <- indata.chunks[[chunk]]
        
      return(ret)
    }) )
    loadings.y <- do.call( c, lapply( names(loadings.condensed.y), function(chunk)
    {
      ret <- rep(loadings.condensed.y[chunk], length(indata.chunks[[chunk]]) )
      names(ret) <- indata.chunks[[chunk]]
      
      return(ret)
    }) )
  }

  colmeans <- colMeans(indata)
  tick.factor <- seq(-2, 2, length = somSize)
  
  
  weightMatrix <- t( sapply( 1:somSize^2, function(i)
  {
    xi <- (i - 1)%%somSize + 1
    yi <- (i - 1)%/%somSize + 1
    return( colmeans + tick.factor[xi] * loadings.x + tick.factor[yi] * loadings.y )
  } ) )
  
  return(weightMatrix)  
}


som.linear.init <- function(indata, somSize, batch)
{
  if( missing("batch") || length(batch) != ncol(indata) ) batch <- rep(1,ncol(indata))
   
  weightMatrix <- do.call( cbind, lapply( unique(batch), function(batch.x)
  {
    indata.sub <- indata[,which(batch==batch.x)]
    indata.sub.noNA <- indata.sub[ apply(indata.sub,1, function(x){all(!is.na(x))}) , ]
    
    som.linear.init.subdata(indata.sub.noNA, somSize )
  }) )
  if( any(batch!=1) ) weightMatrix <- weightMatrix[ , colnames(indata) ]  

  return(weightMatrix)
}