pipeline.sampleSimilarityAnalysisSOM <- function()
{

  filename <- file.path(paste(files.name, "- Results"), "Sample Similarity Analysis", "Sample SOM.pdf")
  util.info("Writing:", filename)
  pdf(filename, 21/2.54, 21/2.54, useDingbats=FALSE)

  for (i in 1:2 )
  {
    d <- if( i == 1 ) get(paste("spot.list.",preferences$standard.spot.modules,sep=""))$spotdata else metadata
    n <- if( i == 1 ) paste( "module data (", preferences$standard.spot.modules, ")") else "metagene data"
    
    secLvlSom <- som(t(d), xdim=preferences$dim.2ndLvlSom, ydim=preferences$dim.2ndLvlSom)
    
    ##### Plot SampleSOM #####
    
    par(mar=c(1, 1, 1, 1))
    xl <- c(min(secLvlSom$visual[,"x"])-1.2, max(secLvlSom$visual[,"x"])+1.2)
    yl <- c(min(secLvlSom$visual[,"y"])-1.2, max(secLvlSom$visual[,"y"])+1.2)
  
    plot( 0, type="n", axes=FALSE, xlab="", ylab="", xlim=xl, ylim=yl, xaxs="i", yaxs="i",main=paste("Sample SOM on",n), cex.main=1)
  
    if (length(unique(group.labels)) > 1)
    {
      legend("topright", as.character(unique(group.labels)), cex=0.5,
             text.col=groupwise.group.colors, bg="white")
    }
  
    for (j in 1:nrow(secLvlSom$code.sum))
    {
  
      which.samples <-
        intersect(which(secLvlSom$visual[,"x"] == secLvlSom$code.sum[j,"x"]),
                  which(secLvlSom$visual[,"y"] == secLvlSom$code.sum[j,"y"]))
  
      if (!is.na(which.samples[1]))
      {
  
        which.samples <- which.samples[1:min(9, length(which.samples))]
  
        x.seq <- c(0,0.3,0,-0.3,0,-0.3,0.3,-0.3,0.3)[seq_along(which.samples)]
        y.seq <- c(0,0,0.3,0,-0.3,-0.3,-0.3,0.3,0.3)[seq_along(which.samples)]
  
        points(secLvlSom$visual[which.samples[1], "x"]+x.seq,
               secLvlSom$visual[which.samples[1], "y"]+y.seq,
               pch=16, col=group.colors[which.samples], cex=1.5)
  
        points(secLvlSom$visual[which.samples[1], "x"]+x.seq,
               secLvlSom$visual[which.samples[1], "y"]+y.seq,
               pch=1, col="gray20", cex=1.5, lwd=1)
      }
    }
    box()
  
    
    if( ncol(d) < 1000 )
    {
      plot( 0, type="n", axes=FALSE, xlab="", ylab="", xlim=xl, ylim=yl, xaxs="i", yaxs="i",main=paste("Sample SOM on",n), cex.main=1)
    
      if (length(unique(group.labels)) > 1)
      {
        legend("topright", as.character(unique(group.labels)), cex=0.5,
               text.col=groupwise.group.colors, bg="white")
      }
    
      for (j in 1:nrow(secLvlSom$code.sum))
      {
    
        which.samples <-
          intersect(which(secLvlSom$visual[,"x"] == secLvlSom$code.sum[j,"x"]),
                    which(secLvlSom$visual[,"y"] == secLvlSom$code.sum[j,"y"]))
    
        if (!is.na(which.samples[1]))
        {
    
          which.samples <- which.samples[1:min(9, length(which.samples))]
    
          x.seq <- c(0,0.3,0,-0.3,0,-0.3,0.3,-0.3,0.3)[seq_along(which.samples)]
          y.seq <- c(0,0,0.3,0,-0.3,-0.3,-0.3,0.3,0.3)[seq_along(which.samples)]
    
          points(secLvlSom$visual[which.samples[1], "x"]+x.seq,
                 secLvlSom$visual[which.samples[1], "y"]+y.seq,
                 pch=16, col=group.colors[which.samples], cex=1.5)
    
          points(secLvlSom$visual[which.samples[1], "x"]+x.seq,
                 secLvlSom$visual[which.samples[1], "y"]+y.seq,
                 pch=1, col="gray20", cex=1.5, lwd=1)
    
          text(secLvlSom$visual[which.samples[1], "x"]+x.seq,
               secLvlSom$visual[which.samples[1], "y"]+y.seq,
               colnames(indata)[which.samples], col="gray20", cex=0.6)
        }
      }
    
      box()
    }
    
    
    ##### Plot SampleSOM with expression portraits ######
    par(mar=c(1,1,1,1))
    xl <- c(0,21)
    yl <- c(0,21)
  
    plot( 0, type="n", axes=FALSE, xlab="", ylab="", xlim=xl, ylim=yl, xaxs="i", yaxs="i",main=paste("Sample SOM on",n), cex.main=1)  
    
    if (length(unique(group.labels)) > 1)
    {
      legend("topright", as.character(unique(group.labels)), cex=0.5,
             text.col=groupwise.group.colors, bg="white")
    }
  
    coord.bins <- as.numeric( cut( c(1:preferences$dim.2ndLvlSom), breaks=20 ) )
    for (j in 1:nrow(secLvlSom$code.sum))
    {
      
      which.samples <-
        intersect(which(coord.bins[secLvlSom$visual[,"x"]+1] == coord.bins[secLvlSom$code.sum[j,"x"]+1]),
                  which(coord.bins[secLvlSom$visual[,"y"]+1] == coord.bins[secLvlSom$code.sum[j,"y"]+1]))
  
      if (!is.na(which.samples[1]))
      {
        m <- matrix(metadata[, which.samples[1]], preferences$dim.1stLvlSom, preferences$dim.1stLvlSom)
  
        if (max(m) - min(m) != 0)
        {
          m <- (m - min(m)) / (max(m) - min(m)) * 999
        }
  
        m <- cbind(apply(m, 1, function(x){x}))[nrow(m):1,]
        x <- pixmapIndexed(m , col=color.palette.portraits(1000), cellres=10)
  
        addlogo(x,
                coord.bins[secLvlSom$visual[which.samples[1],"x"]+1]+c(-0.45,0.455),
                coord.bins[secLvlSom$visual[which.samples[1],"y"]+1]+c(-0.45,0.45))
  
        which.samples <- which.samples[1:min(9, length(which.samples))]
  
        x.seq <- c(0,0.3,0,-0.3,0,-0.3,0.3,-0.3,0.3)[seq_along(which.samples)]
        y.seq <- c(0,0,0.3,0,-0.3,-0.3,-0.3,0.3,0.3)[seq_along(which.samples)]
  
        points(coord.bins[secLvlSom$visual[which.samples[1],"x"]+1]+x.seq,
               coord.bins[secLvlSom$visual[which.samples[1],"y"]+1]+y.seq,
               pch=16, col=group.colors[which.samples], cex=1.2)
  
        points(coord.bins[secLvlSom$visual[which.samples[1],"x"]+1]+x.seq,
               coord.bins[secLvlSom$visual[which.samples[1],"y"]+1]+y.seq,
               pch=1, col="gray20", cex=1.2, lwd=1)
      }
    }
  
    box()
    
  }  
    
  dev.off()
}
