pipeline.patAssignment <- function()
{
  spot.list <- get(paste("spot.list.",preferences$standard.spot.modules,sep=""))

  thresh.global <- sd(as.vector(spot.list$spotdata))
  
  pat.labels <<- apply( spot.list$spotdata > thresh.global, 2, function(x)
  {
    paste( names(x)[which(x)], collapse=" " )
  } )
  
  # join small pats into their precursors
  
  spot.counts <- rowSums( spot.list$spotdata > thresh.global )
  skip.pats <- ""
  while( sort(table(pat.labels[which(!pat.labels%in%skip.pats)]))[1] < length(pat.labels)*0.01 )
  {
    pat.to.merge <- names(sort(table(pat.labels[which(!pat.labels%in%skip.pats)]))[1])
    if( length(strsplit(pat.to.merge," ")[[1]] > 1) )
    {
      least.freq.spot <- names(sort(spot.counts[ strsplit(pat.to.merge," ")[[1]] ])[1])
      pat.after.merge <- sub(least.freq.spot,"",pat.to.merge)
      pat.after.merge <- sub("  "," ",pat.after.merge)
      pat.after.merge <- sub("^ | $","",pat.after.merge)  
      pat.labels[which(pat.labels==pat.to.merge)] <<- pat.after.merge
    } else
    {
      skip.pats <- c(skip.pats,pat.to.merge)
    }
  } 
  
  pat.labels[which(pat.labels=="")] <<- "none"
}