modules.profiles <- function(spot.list, main, path)
{
  sd.theshold = sd(spot.list$spotdata)
  
  pdf(path, 29.7/2.54, 21/2.54)
  
  #### spot profiles ####
  layout(matrix(1:15, 5, 3, byrow=TRUE), widths=c(1,4,1))
  
  for (i in 1:length(spot.list$spots))
  {
    if (main %in% c("Sample-Underexpression"))
    {
      samples <- which(spot.list$spotdata[i,] < -sd.theshold)
    } else
    {
      samples <- which(spot.list$spotdata[i,] > sd.theshold)
    }
    
    par(mar=c(0.5,3,0.5,1))
    mask <- spot.list$spots[[i]]$mask
    image(matrix(mask, preferences$dim.1stLvlSom, preferences$dim.1stLvlSom),
          axes=FALSE, col ="darkgreen")
      axis(2, 0.95, names(spot.list$spots[i]), las=2, tick=FALSE, cex.axis=1.6)
      box()
    
    par(mar=c(0.5,3,0.5,1))
    
    barplot.x <- barplot(spot.list$spotdata[i,], names.arg=NA,
                         col=group.colors, main="", las=2, cex.main=1,
                         cex.lab=1, cex.axis=1, cex.names=0.8,
                         border=if (ncol(indata) < 80) "black" else NA,
                         ylim=range(spot.list$spotdata)*c(1,1.4))
    
    abline(h=sd.theshold*c(-1,1), lty=2, col="gray")

    label.string <- rep(".",length(unique(group.labels)))
    label.string[  which( tapply( spot.list$spotdata[i,], group.labels, mean )[unique(group.labels)] > sd.theshold )  ] = "+"
    label.string[  which( tapply( spot.list$spotdata[i,], group.labels, mean )[unique(group.labels)] < -sd.theshold )  ] = "-"
    label.x <- tapply(barplot.x, group.labels, mean)[unique(group.labels)]
    text(label.x, max(spot.list$spotdata)*1.2, label.string, col=groupwise.group.colors,cex=2.5)
    
    par(mar=c(0.5,0,0.5,0))
    plot(0, type="n", axes=FALSE, xlab="", ylab="", xlim=c(0,1), ylim=c(0,1), xaxs="i", yaxs="i")
    text(0, 0.9, paste("# samples :", length(samples)), adj=0)
    
    if (length(samples) > 0)
    {
      t <- table(as.vector(group.labels[samples]))[unique(group.labels)]
      t <- t[!is.na(t)]
      
      for (ii in seq_along(t))
      {
        t[ii] <- t[ii] / table(as.vector(group.labels))[names(t[ii])]
      }
      
      t <- round(100 * t)
      
      t2 <- table(as.vector(group.labels[samples]))[unique(group.labels)]
      t2 <- t2[!is.na(t2)]
      
      for (ii in seq_along(t))
      {
        text(0, 0.8  - ii*0.1,
             paste(names(t)[ii], ": # =", t2[ii], " -> ", t[ii], "%"),
             adj=0, cex=1, col=group.colors[which(group.labels == names(t)[ii])[1]])
      }
    }
  }
  
  ### plot profile boxplots
  layout(matrix(1:8, 4, 2, byrow=TRUE), widths=c(1,4))
  ylim <- range( unlist( by( t(spot.list$spotdata), group.labels, range )[unique(group.labels)] ) )
  
  for (i in 1:length(spot.list$spots))
  {
    par(mar=c(2,3,0.5,1))
    
    image(matrix(spot.list$spots[[i]]$mask, preferences$dim.1stLvlSom, preferences$dim.1stLvlSom),
          axes=FALSE, col ="darkgreen")
    
    axis(2, 0.95, names(spot.list$spots[i]), las=2, tick=FALSE, cex.axis=1.6)
    box()
    
    par(mar=c(2,3,0.5,1))
    
    boxplot(  tapply( spot.list$spotdata[i,], group.labels, c )[unique(group.labels)],
              col=groupwise.group.colors, las=1, ylim=ylim )
    
    abline(h=sd.theshold*c(-1,0,1), lty=2, col="gray")
    
    label.string <- rep(".",length(unique(group.labels)))
    label.string[  which( tapply( spot.list$spotdata[i,], group.labels, mean )[unique(group.labels)] > sd.theshold )  ] = "+"
    label.string[  which( tapply( spot.list$spotdata[i,], group.labels, mean )[unique(group.labels)] < -sd.theshold )  ] = "-"
    
    text( 1:length(unique(group.labels)), par("usr")[4]-(par("usr")[4]-par("usr")[3])*0.1, label.string, col=groupwise.group.colors,cex=3)
  }
  
  dev.off()
}