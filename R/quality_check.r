pipeline.qualityCheck <- function(env)
{
  dir.create("Data Overview", showWarnings=FALSE)

  plot.poly.density = function(subset=1:ncol(env$indata), col="#BEBEBE", main="", add=FALSE)
  {
    col.t=paste(col,"50",sep="")

    if (add) par(new=TRUE)

    plot(densities.y[subset[1],], type="n", xlim=range(densities.x), ylim=c(0, max(densities.y)), col="white", main=main, xlab="log Expression", ylab="Density")

    lines(densities.x, colMeans(densities.y[subset,,drop=FALSE]), lwd=2, col=col)

    polygon(
      x= c(densities.x, rev(densities.x)),
      y= c(colMeans(densities.y[subset,,drop=FALSE])+ apply(densities.y[subset,,drop=FALSE], 2, sd), rev(colMeans(densities.y[subset,,drop=FALSE]) - apply(densities.y[subset,,drop=FALSE], 2, sd))),
      col=col.t, border=col)
  }

  mi = min(env$indata); ma = max(env$indata)
  densities = apply(env$indata,2,function(x) density(x, from=mi, to=ma))
  densities.y = t(sapply(densities, function(x) x$y))
  densities.x = densities[[1]]$x

  filename <- file.path("Data Overview", "Data Distribution.pdf")
  util.info("Writing:", filename)
  pdf(filename, 29.7/2.54, 21/2.54, useDingbats=FALSE)

  par(mfrow=c(1,2))
  plot(densities.x, densities.y[1,], main="Input data distribution", xlim=range(env$indata), ylim=range(densities.y), type="l", xlab="log Expression", ylab="Density")
    dummy=sapply(densities, lines)

  for (i in seq_along(unique(env$group.labels)))
    plot.poly.density(subset=which(env$group.labels==unique(env$group.labels)[i]), col=env$groupwise.group.colors[i], main="", add=ifelse(i==1,FALSE,TRUE))

  par(mfrow=c(2,3))
  for (i in seq_along(unique(env$group.labels)))
    plot.poly.density(subset=which(env$group.labels==unique(env$group.labels)[i]), col=env$groupwise.group.colors[i], main=unique(env$group.labels)[i])

	
	indata.sample.sd <- apply(env$indata, 2, sd)
	indata.sample.mean <- colMeans(env$indata)

	Q13 <- quantile( indata.sample.mean, c(0.25,0.75) )
	IQR1.mean <- c( Q13[1] - 1*diff( Q13 ), Q13[2] + 1*diff( Q13 ) )
	IQR3.mean <- c( Q13[1] - 3*diff( Q13 ), Q13[2] + 3*diff( Q13 ) )
	Q13 <- quantile( indata.sample.sd, c(0.25,0.75) )
	IQR1.sd <- c( Q13[1] - 1*diff( Q13 ), Q13[2] + 1*diff( Q13 ) )
	IQR3.sd <- c( Q13[1] - 3*diff( Q13 ), Q13[2] + 3*diff( Q13 ) )

	outlier <- names( which( indata.sample.mean < IQR1.mean[1] |
														 indata.sample.mean > IQR1.mean[2] |
														 indata.sample.sd < IQR1.sd[1] |
														 indata.sample.sd > IQR1.sd[2] ) )

	par(mfrow=c(1,1), mar=c(5,4,3,2))
	plot( indata.sample.mean, indata.sample.sd, pch=16, col=env$group.colors, xlab="mean expression level", ylab="standard deviation of expression" )
		abline( v=IQR1.mean, col="gray20", lty=2 )
		abline( v=IQR3.mean, col="gray20", lty=3 )
		abline( h=IQR1.sd, col="gray20", lty=2 )
		abline( h=IQR3.sd, col="gray20", lty=3 )
		legend("topright",c("1x IQR","3x IQR"),lty=c(2,3),col="gray20")
	if(length(outlier)>0) text(indata.sample.mean[outlier],indata.sample.sd[outlier]+diff(range(indata.sample.sd))*0.017,outlier)
		
		
  par(mfrow=c(2,1),mar=c(5,3,3,2))

  barplot(indata.sample.mean, col=env$group.colors, main="Sample mean expression", names.arg=if (ncol(env$indata)<80) colnames(env$indata) else rep("",ncol(env$indata)), las=2, cex.main=1, cex.lab=1, cex.axis=0.8, cex.names=0.6, border = ifelse(ncol(env$indata)<80,"black",NA), ylim=range(indata.sample.mean)*c(0.99,1.01), xpd=FALSE)
    box()
	

  if (length(unique(env$group.labels)) > 1)
  {
    mean.boxes = by(indata.sample.mean, env$group.labels, c)[unique(env$group.labels)]
    par(mar=c(5,3,0,2))
    boxplot(mean.boxes, col=env$groupwise.group.colors, las=2, main="", cex.main=1, cex.axis=0.8, xaxt="n")
    axis(1, seq_along(env$groupwise.group.colors), unique(env$group.labels), las=2, cex.axis=0.8)
  }


  par(mar=c(5,3,3,2))

  barplot(indata.sample.sd, col=env$group.colors, main="Sample expression standard deviation", names.arg=if (ncol(env$indata)<80) colnames(env$indata) else rep("",ncol(env$indata)), las=2, cex.main=1, cex.lab=1, cex.axis=0.8, cex.names=0.6, border = ifelse(ncol(env$indata)<80,"black",NA), ylim=range(indata.sample.sd)*c(0.99,1.01), xpd=FALSE)
  box()

  if (length(unique(env$group.labels)) > 1)
  {
    mean.boxes = by(indata.sample.sd, env$group.labels, c)[unique(env$group.labels)]
    par(mar=c(5,3,0,2))
    boxplot(mean.boxes, col=env$groupwise.group.colors, las=2, main="", cex.main=1, cex.axis=0.8, xaxt="n")
    axis(1, seq_along(env$groupwise.group.colors), unique(env$group.labels), las=2, cex.axis=0.8)
  }
  
  
  pipeline.affymetrixQualityCheck(env)


  dev.off()
}
