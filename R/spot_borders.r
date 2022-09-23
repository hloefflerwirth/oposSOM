spot.borders <- function(spot.list,...)
{
  spot.mask <- matrix(spot.list$overview.mask,sqrt(length(spot.list$overview.mask)) )
  spot.mask[which(is.na(spot.mask))] <- 0
  
  par(new=T)
  # plot(0, type="n", xlim=c(1,env$preferences$dim), ylim=par("usr")[1:2], xlab="",ylab="", axes=F, xaxs="i", yaxs="i"  )
  plot(0, type="n", xlim=c(.5,env$preferences$dim.1stLvlSom+.5), ylim=c(.5,env$preferences$dim.1stLvlSom+.5), xlab="",ylab="", axes=F, xaxs="i", yaxs="i" )
  
  for( x in 1:nrow(spot.mask) )
    for( y in 1:nrow(spot.mask) )
    {
      nei.r <- oposSOM::get.neighbors(x,y,nrow(spot.mask))$r
      nei.d <- oposSOM::get.neighbors(x,y,nrow(spot.mask))$d
      
      if( !is.null(nei.r) && spot.mask[x,y] != spot.mask[nei.r[1],nei.r[2]] ) lines(c(x+0.5,x+0.5),c(y-0.5,y+0.5),...)
      if( !is.null(nei.d) && spot.mask[x,y] != spot.mask[nei.d[1],nei.d[2]] ) lines(c(x-0.5,x+0.5),c(y-0.5,y-0.5),...)
    }
}