pipeline.groupAssignment <- function()
{
  if (length(unique(group.labels)) < 2)
  {
    group.bootstrap.score <<- rep(100, ncol(indata))
    return()
  }


  PCM <- cor( metadata )
  diag(PCM) <- NA
  
  group.silhouette.coef <<- sapply( seq(ncol(metadata)), function(i)
  {
    mean.group.correlations <- tapply( PCM[,i], group.labels, mean, na.rm=TRUE )
    
    cor.m.A <- mean.group.correlations[ group.labels[i] ]
    cor.m.B <- max( mean.group.correlations[ which( names(mean.group.correlations) != group.labels[i] ) ] )
    
    return( cor.m.A - cor.m.B )
  } )
  names(group.silhouette.coef) <<- colnames(metadata)

  group.silhouette.coef[which(is.nan(group.silhouette.coef))] <<- 0

}
