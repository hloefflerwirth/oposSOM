

psf.overview.heatmaps <- function(psf.results, output.path, group.colors, color.palette )
{
  
  filename <- file.path(output.path, "0verview Heatmaps.pdf")
  util.info("Writing:", filename)
  pdf(filename, 21/2.54, 21/2.54)
  
  
  mean.psf.matrix <- t( sapply( psf.results, function(x) colMeans( x ) ) )
  mean.psf.matrix <- mean.psf.matrix[ order(apply(mean.psf.matrix,1,var),decreasing=TRUE)[1:min(20,nrow(mean.psf.matrix))] , ,drop=F]
  
  heatmap(x=log1p(log1p(mean.psf.matrix)), cex.main=2,
               col=color.palette(1000),scale="r",
               mar=c(10,20), ColSideColors=group.colors, cexDend=0.6 )
  
  par(new=TRUE, mar=c(3.8,29.5,36.2,2.5))
  
  image(matrix(c(1:1000), 1000, 1), axes=FALSE,
        col=color.palette(1000))
  
  box()
  axis(1, c(0,1), round( range( log10(mean.psf.matrix)) ), cex.axis=1)
  mtext(bquote("log"[10] ~ "s"), side=1, cex=1, line=1)
  mtext("mean node signal", cex=1.4, line=33)
  
  
  heatmap(x=log1p(log1p(mean.psf.matrix)), cex.main=2,
               col=color.palette(1000),scale="r",
               mar=c(10,20), ColSideColors=group.colors, Colv=NA, cexDend=0.6)
  
  par(new=TRUE, mar=c(3.8,29.5,36.2,2.5))
  
  image(matrix(c(1:1000), 1000, 1), axes=FALSE,
        col=color.palette(1000))
  
  box()
  axis(1, c(0,1), round( range( log10(mean.psf.matrix)) ), cex.axis=1)
  mtext(bquote("log"[10] ~ "s"), side=1, cex=1, line=1)
  mtext("mean node signal", cex=1.4, line=33)
  
  
  
  # max.psf.matrix <- t( sapply( psf.results, function(x) sapply( x, function(y) if(length(y$signal.at.sinks)>0) max(y$signal.at.sinks) else 1 )  ) )
  # max.psf.matrix <- max.psf.matrix[ order(apply(max.psf.matrix,1,var),decreasing=TRUE)[1:(nrow(max.psf.matrix)/2)] , ]
  # 
  # heatmap(x=log1p(log1p(max.psf.matrix)), cex.main=2,
  #            col=colorRampPalette(c("blue4","blue","gray90","orange","red4"))(1000),scale="r",
  #            mar=c(10,20), ColSideColors=group.colors, cexDend=0.6 )
  # 
  # par(new=TRUE, mar=c(3.8,29.5,36.2,2.5))
  # 
  # image(matrix(c(1:1000), 1000, 1), axes=FALSE,
  #       col=colorRampPalette(c("blue4","blue","gray90","orange","red4"))(1000))
  # 
  # box()
  # axis(1, c(0,1), round( range( log10(max.psf.matrix)) ), cex.axis=1)
  # mtext(bquote("log"[10] ~ "s"), side=1, cex=1, line=1)
  # mtext("signal maximum", cex=1.4, line=33)
  # 
  # 
  # heatmap(x=log1p(log1p(max.psf.matrix)), cex.main=2,
  #            col=colorRampPalette(c("blue4","blue","gray90","orange","red4"))(1000),scale="r",
  #            mar=c(10,20), ColSideColors=group.colors, Colv=NA, cexDend=0.6)
  # 
  # par(new=TRUE, mar=c(3.8,29.5,36.2,2.5))
  # 
  # image(matrix(c(1:1000), 1000, 1), axes=FALSE,
  #       col=colorRampPalette(c("blue4","blue","gray90","orange","red4"))(1000))
  # 
  # box()
  # axis(1, c(0,1), round( range( log10(max.psf.matrix)) ), cex.axis=1)
  # mtext(bquote("log"[10] ~ "s"), side=1, cex=1, line=1)
  # mtext("signal maximum", cex=1.4, line=33)
  
  
  dev.off() 
  
}














