pipeline.qualityCheck <- function()
{
  dir.create(output.paths["LPE"], showWarnings=FALSE)

  plot.poly.density = function(subset=1:ncol(indata), col="#BEBEBE", main="", add=FALSE)
  {
    col.t=paste(col,"50",sep="")

    if (add) par(new=TRUE)

    plot(densities.y[subset[1],], type="n", xlim=range(densities.x), ylim=c(0, max(densities.y)), col="white", main=main, xlab="log Expression", ylab="Density")

    lines(densities.x, colMeans(densities.y[subset,,drop=FALSE]), lwd=2, col=col)

    polygon(
      x= c(densities.x, rev(densities.x)),
      y= c(colMeans(densities.y[subset,,drop=FALSE])+ apply(densities.y[subset,,drop=FALSE], 2, sd), rev(colMeans(densities.y[subset,,drop=FALSE]) - apply(densities.y[subset,,drop=FALSE], 2, sd))),
      col=col.t, border=col
)
  }

  mi = min(indata); ma = max(indata)
  densities = apply(indata,2,function(x) density(x, from=mi, to=ma))
  densities.y = t(sapply(densities, function(x) x$y))
  densities.x = densities[[1]]$x

  pdf(paste(output.paths["LPE"],"/data distribution.pdf", sep=""), 29.7/2.54, 21/2.54)

  par(mfrow=c(1,2))

  plot(densities.x, densities.y[1,], main="Input data distribution", xlim=range(indata), ylim=range(densities.y), type="l", xlab="log Expression", ylab="Density")
    dummy=sapply(densities, lines)

  for (i in 1:length(unique(group.labels)))
    plot.poly.density(subset=which(group.labels==unique(group.labels)[i]), col=groupwise.group.colors[i], main="", add=ifelse(i==1,FALSE,TRUE))

  par(mfrow=c(2,3))
  for (i in 1:length(unique(group.labels)))
    plot.poly.density(subset=which(group.labels==unique(group.labels)[i]), col=groupwise.group.colors[i], main=unique(group.labels)[i])

  par(mfrow=c(2,1))
  par(mar=c(5,3,3,2))

  barplot(indata.sample.mean, col=group.colors, main="Sample mean expression", names.arg=if (ncol(indata)<80) colnames(indata) else rep("",ncol(indata)), las=2, cex.main=1, cex.lab=1, cex.axis=0.8, cex.names=0.6, border = ifelse(ncol(indata)<80,"black",NA), ylim=range(indata.sample.mean)*c(0.99,1.01), xpd=FALSE)
    box()


  if (length(unique(group.labels)) > 1)
  {
    mean.boxes = by(indata.sample.mean, group.labels, c)[unique(group.labels)]
    par(mar=c(5,3,0,2))
    boxplot(mean.boxes, col=groupwise.group.colors, las=2, main="", cex.main=1, cex.axis=0.8, xaxt="n")
    axis(1, 1:length(groupwise.group.colors), unique(group.labels), las=2, cex.axis=0.8)
  }



  indata.sample.var <- apply(indata, 2, var)

  par(mar=c(5,3,3,2))

  barplot(indata.sample.var, col=group.colors, main="Sample expression variance", names.arg=if (ncol(indata)<80) colnames(indata) else rep("",ncol(indata)), las=2, cex.main=1, cex.lab=1, cex.axis=0.8, cex.names=0.6, border = ifelse(ncol(indata)<80,"black",NA), ylim=range(indata.sample.var)*c(0.99,1.01), xpd=FALSE)
  box()


  if (length(unique(group.labels)) > 1)
  {
    mean.boxes = by(indata.sample.var, group.labels, c)[unique(group.labels)]
    par(mar=c(5,3,0,2))
    boxplot(mean.boxes, col=groupwise.group.colors, las=2, main="", cex.main=1, cex.axis=0.8, xaxt="n")
    axis(1, 1:length(groupwise.group.colors), unique(group.labels), las=2, cex.axis=0.8)
  }


  dev.off()
}
