pipeline.patAssignment <- function(env)
{
  find.next.merge.pat <- function( pat.labels, skip.pats, spot.counts )
  {
    tab <- table(pat.labels[which(!pat.labels%in%skip.pats)])
    
    candidates <- names( which( tab == min(tab) ) )
    candidates <- candidates[ which( nchar(candidates) == max(nchar(candidates)) ) ]
    
    candidates.counts <- sapply( candidates, function(x) sum( spot.counts[ strsplit(x," ")[[1]] ] )    )
    
    return( names(which.min(candidates.counts)) )    
  }
	
	spot.list <- env[[paste("spot.list.", env$preferences$standard.spot.modules,sep="")]]

  thresh.global <- sd(as.vector(spot.list$spotdata))
  
  env$pat.labels <- apply( spot.list$spotdata > thresh.global, 2, function(x)
  {
    paste( names(x)[which(x)], collapse=" " )
  } )
  
  # join small pats into their precursors
  
  spot.counts <- rowSums( spot.list$spotdata > thresh.global )
  skip.pats <- ""
  if(any(env$pat.labels!=""))
    while( sort(table(env$pat.labels[which(!env$pat.labels%in%skip.pats)]))[1] < length(env$pat.labels)*0.01 )
    {
      pat.to.merge <- find.next.merge.pat( env$pat.labels, skip.pats, spot.counts )
#      if( length(strsplit(pat.to.merge," ")[[1]] ) > 1 ) 
#      {
        least.freq.spot <- names(sort(spot.counts[ strsplit(pat.to.merge," ")[[1]] ])[1])
        pat.after.merge <- sub(least.freq.spot,"",pat.to.merge)
        pat.after.merge <- sub("  "," ",pat.after.merge)
        pat.after.merge <- sub("^ | $","",pat.after.merge)  
        env$pat.labels[which(env$pat.labels==pat.to.merge)] <- pat.after.merge
 #     } else
 #     {
 #       skip.pats <- c(skip.pats,pat.to.merge)
 #     }
    } 
  
  env$pat.labels[which(env$pat.labels=="")] <- "none"
  return(env)
}