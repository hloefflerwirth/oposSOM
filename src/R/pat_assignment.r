pipeline.patAssignment <- function()
{
  find.next.merge.pat <- function( pat.labels, skip.pats, spot.counts )
  {
    tab <- table(pat.labels[which(!pat.labels%in%skip.pats)])
    
    candidates <- names( which( tab == min(tab) ) )
    candidates <- candidates[ which( nchar(candidates) == max(nchar(candidates)) ) ]
    
    candidates.counts <- sapply( candidates, function(x) sum( spot.counts[ strsplit(x," ")[[1]] ] )    )
    
    return( names(which.min(candidates.counts)) )    
  }
	
	
  spot.list <- get(paste("spot.list.",preferences$standard.spot.modules,sep=""))

  thresh.global <- sd(as.vector(spot.list$spotdata))
  
  pat.labels <<- apply( spot.list$spotdata > thresh.global, 2, function(x)
  {
    paste( names(x)[which(x)], collapse=" " )
  } )
  
  # join small pats into their precursors
  
  spot.counts <- rowSums( spot.list$spotdata > thresh.global )
  skip.pats <- ""
  if(any(pat.labels!=""))
    while( sort(table(pat.labels[which(!pat.labels%in%skip.pats)]))[1] < length(pat.labels)*0.01 )
    {
      pat.to.merge <- find.next.merge.pat( pat.labels, skip.pats, spot.counts )
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