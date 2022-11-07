
color.palette.discrete <- function(n)   # color scheme from https://personal.sron.nl/~pault/colourschemes.pdf
{
  color.set <- c( "#E8ECFB", "#D9CCE3", "#D1BBD7", "#CAACCB", "#BA8DB4",
                  "#AE76A3", "#AA6F9E", "#994F88", "#882E72", "#1965B0",
                  "#437DBF", "#5289C7", "#6195CF", "#7BAFDE", "#4EB265",
                  "#90C987", "#CAE0AB", "#F7F056", "#F7CB45", "#F6C141",
                  "#F4A736", "#F1932D", "#EE8026", "#E8601C", "#E65518",
                  "#DC050C", "#A5170E", "#72190E", "#42150A" )
  
  n <- ceiling(n)
    
  if (n == 1) return( color.set[c(10)] )
  if (n == 2) return( color.set[c(10,26)] )
  if (n == 3) return( color.set[c(10,18,26)] )
  if (n == 4) return( color.set[c(10,15,18,26)] )
  if (n == 5) return( color.set[c(10,14,15,18,26)] )
  
  if (n == 6) return( color.set[c(10,14,15,17,18,26)] )
  if (n == 7) return( color.set[c(9,10,14,15,17,18,26)] )
  if (n == 8) return( color.set[c(9,10,14,15,17,18,23,26)] )
  if (n == 9) return( color.set[c(9,10,14,15,17,18,23,26,28)] )
  if (n == 10) return( color.set[c(9,10,14,15,17,18,21,24,26,28)] )
                                 
  if (n == 11) return( color.set[c(9,10,12,14,15,17,18,21,24,26,28)] )
  if (n == 12) return( color.set[c(3,6,9,10,12,14,15,17,18,21,24,26)] )
  if (n == 13) return( color.set[c(3,6,9,10,12,14,15,16,17,18,21,24,26)] )
  if (n == 14) return( color.set[c(3,6,9,10,12,14,15,16,17,18,20,22,24,26)] )
  if (n == 15) return( color.set[c(3,6,9,10,12,14,15,16,17,18,20,22,24,26,28)] )
  
  if (n == 16) return( color.set[c(3,5,7,9,10,12,14,15,16,17,18,20,22,24,26,28)] )
  if (n == 17) return( color.set[c(3,5,7,8,9,10,12,14,15,16,17,18,20,22,24,26,28)] )
  if (n == 18) return( color.set[c(3,5,7,8,9,10,12,14,15,16,17,18,20,22,24,26,27,28)] )
  if (n == 19) return( color.set[c(2,4,5,7,8,9,10,12,14,15,16,17,18,20,22,24,26,27,28)] )
  if (n == 20) return( color.set[c(2,4,5,7,8,9,10,11,13,14,15,16,17,18,20,22,24,26,27,28)] )
  
  if (n == 21) return( color.set[c(2,4,5,7,8,9,10,11,13,14,15,16,17,18,19,21,22,24,26,27,28)] )
  if (n == 22) return( color.set[c(2,4,5,7,8,9,10,11,13,14,15,16,17,18,19,21,22,24,26,27,28,29)] )
  # if (n == 23) return( color.set[c(1,2,4,5,7,8,9,10,11,13,14,15,16,17,18,19,21,22,24,26,27,28,29)] )

  if (n >= 23) return( colorRampPalette(color.set[-1])(n) )
}

color.palette.interlace <- function(n)
{
  col <- color.palette.discrete(n)
  o <- do.call( c, lapply( 1:floor(n/2), function(i) c(i,i+ceiling(n/2)) ) )
  if( length(o)==n-1 ) o <- c( o, ceiling(n/2) )
  return(col[o])
}

color.palette.portraits <- colorRampPalette(c("darkblue","blue","lightblue3","green3","yellow2","red2","darkred"))

# color.palette.heatmaps <- colorRampPalette(c("#2C7BB6","#64A4CC","#9CCEE3","#C6E5DB","#ECF6C8","#FEEDAA","#FDC980","#F89D59","#E75B3A","#D7191C"))
color.palette.heatmaps <- colorRampPalette(c("#364B9A","#4A7BB7","#6EA6CD","#98CAE1","#C2E4EF","#EAECCC","#FEDA8B","#FDB366","#F67E4B","#DD3D2D","#A50026"))



