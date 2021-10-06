pipeline.patAssignment <- function(env)
{
  find.next.merge.pat <- function( pat.labels, spot.counts )
  {
    tab <- table(pat.labels)
    
    candidates <- names( which( tab == min(tab) ) )
    candidates <- candidates[ which( nchar(candidates) == max(nchar(candidates)) ) ]
    
    candidates.counts <- sapply( candidates, function(x) sum( spot.counts[ strsplit(x," ")[[1]] ] )    )
    
    return( names(which.min(candidates.counts)) )    
  }
	
	spot.list <- env[[paste("spot.list.", env$preferences$standard.spot.modules,sep="")]]

  thresh.global <- sd(as.vector(spot.list$spotdata))
  spot.counts <- rowSums( spot.list$spotdata > thresh.global )
  spot.counts <- spot.counts[which(spot.counts>0)]
  spot.order <- order(spot.counts,decreasing = T)
  
  env$pat.labels <- apply( spot.list$spotdata, 2, function(x)
  {
    x <- x[ spot.order ] > thresh.global
    return(   paste( names(x)[which(x)], collapse=" " )   )
  } )

  
  # join small pats into their precursors
  
  if(any(env$pat.labels!=""))
    while( sort(table(env$pat.labels))[1] < length(env$pat.labels)*0.01 )
    {
      pat.to.merge <- find.next.merge.pat( env$pat.labels, spot.counts )
      least.freq.spot <- names(sort(spot.counts[ strsplit(pat.to.merge," ")[[1]] ])[1])
      pat.after.merge <- sub(least.freq.spot,"",pat.to.merge)
      pat.after.merge <- sub("  "," ",pat.after.merge)
      pat.after.merge <- sub("^ | $","",pat.after.merge)  
      env$pat.labels[which(env$pat.labels==pat.to.merge)] <- pat.after.merge
    } 
  
  env$pat.labels[which(env$pat.labels=="")] <- "none"
	env$pat.colors <- color.palette.discrete(length(unique(env$pat.labels)))[match(env$pat.labels, unique(env$pat.labels))]
	names(env$pat.colors) <- colnames(env$indata)

  return(env)
}