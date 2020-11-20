pipeline.sampleSimilarityAnalysisSOM <- function(env)
{

  filename <- file.path(paste(env$files.name, "- Results"), "Sample Similarity Analysis", "Sample SOM.pdf")
  util.info("Writing:", filename)
  pdf(filename, 21/2.54, 21/2.54, useDingbats=FALSE)

  for (i in 1:2 )
  {
    d <- if( i == 1 ) env[[paste("spot.list.",env$preferences$standard.spot.modules,sep="")]]$spotdata else env$metadata
    n <- if( i == 1 ) paste( "module data (", env$preferences$standard.spot.modules, ")") else "metagene data"
    
    secLvlSom <- som.linear.init(t(d), env$preferences$dim.2ndLvlSom)
    secLvlSom <- som.training(t(d),secLvlSom)
  
    secLvlSom$feature.coords <- secLvlSom$node.summary[secLvlSom$feature.BMU,c("x","y")]

    ##### Plot SampleSOM #####
    
    par(mar=c(1, 1, 1, 1))
    xl <- c(min(secLvlSom$feature.coords[,"x"])-1.2, max(secLvlSom$feature.coords[,"x"])+1.2)
    yl <- c(min(secLvlSom$feature.coords[,"y"])-1.2, max(secLvlSom$feature.coords[,"y"])+1.2)
  
    plot( 0, type="n", axes=FALSE, xlab="", ylab="", xlim=xl, ylim=yl, xaxs="i", yaxs="i",main=paste("Sample SOM on",n), cex.main=1)
  
    if (length(unique(env$group.labels)) > 1)
    {
      legend("topright", as.character(unique(env$group.labels)), cex=0.5,
             text.col=env$groupwise.group.colors, bg="white")
    }
  
    for (j in 1:nrow(secLvlSom$node.summary))
    {
  
      which.samples <-
        intersect(which(secLvlSom$feature.coords[,"x"] == secLvlSom$node.summary[j,"x"]),
                  which(secLvlSom$feature.coords[,"y"] == secLvlSom$node.summary[j,"y"]))
  
      if (!is.na(which.samples[1]))
      {
  
        which.samples <- which.samples[1:min(9, length(which.samples))]
  
        x.seq <- c(0,0.3,0,-0.3,0,-0.3,0.3,-0.3,0.3)[seq_along(which.samples)]
        y.seq <- c(0,0,0.3,0,-0.3,-0.3,-0.3,0.3,0.3)[seq_along(which.samples)]
  
        points(secLvlSom$feature.coords[which.samples[1], "x"]+x.seq,
               secLvlSom$feature.coords[which.samples[1], "y"]+y.seq,
               pch=16, col=env$group.colors[which.samples], cex=1.5)
  
        points(secLvlSom$feature.coords[which.samples[1], "x"]+x.seq,
               secLvlSom$feature.coords[which.samples[1], "y"]+y.seq,
               pch=1, col="gray20", cex=1.5, lwd=1)
      }
    }
    box()
  
    
    if( ncol(d) < 1000 )
    {
      plot( 0, type="n", axes=FALSE, xlab="", ylab="", xlim=xl, ylim=yl, xaxs="i", yaxs="i",main=paste("Sample SOM on",n), cex.main=1)
    
      if (length(unique(env$group.labels)) > 1)
      {
        legend("topright", as.character(unique(env$group.labels)), cex=0.5,
               text.col=env$groupwise.group.colors, bg="white")
      }
    
      for (j in 1:nrow(secLvlSom$node.summary))
      {
    
        which.samples <-
          intersect(which(secLvlSom$feature.coords[,"x"] == secLvlSom$node.summary[j,"x"]),
                    which(secLvlSom$feature.coords[,"y"] == secLvlSom$node.summary[j,"y"]))
    
        if (!is.na(which.samples[1]))
        {
    
          which.samples <- which.samples[1:min(9, length(which.samples))]
    
          x.seq <- c(0,0.3,0,-0.3,0,-0.3,0.3,-0.3,0.3)[seq_along(which.samples)]
          y.seq <- c(0,0,0.3,0,-0.3,-0.3,-0.3,0.3,0.3)[seq_along(which.samples)]
    
          points(secLvlSom$feature.coords[which.samples[1], "x"]+x.seq,
                 secLvlSom$feature.coords[which.samples[1], "y"]+y.seq,
                 pch=16, col=env$group.colors[which.samples], cex=1.5)
    
          points(secLvlSom$feature.coords[which.samples[1], "x"]+x.seq,
                 secLvlSom$feature.coords[which.samples[1], "y"]+y.seq,
                 pch=1, col="gray20", cex=1.5, lwd=1)
    
          text(secLvlSom$feature.coords[which.samples[1], "x"]+x.seq,
               secLvlSom$feature.coords[which.samples[1], "y"]+y.seq,
               colnames(env$indata)[which.samples], col="gray20", cex=0.6)
        }
      }
    
      box()
    }
    
    
    ##### Polygone-2ndSOM #####
    if (length(unique(env$group.labels)) > 1)
    {
      transparent.group.colors <- sapply(env$groupwise.group.colors, function(x)
      {
        paste(substr(x, 1, 7) , "50", sep="")
      })
      
      names(transparent.group.colors) <- unique(env$group.labels)
      
      plot(secLvlSom$feature.coords[,"x"], -secLvlSom$feature.coords[,"y"],
           type="n", axes=FALSE, xlab="", ylab="", cex=4, col=env$group.colors, pch=16,
           xaxs="i", yaxs="i", xlim=xl, ylim=yl,main="Second level SOM, group outlines", cex.main=1)
      
      for (i in seq_along(unique(env$group.labels)))
      {
        group.member <- which(env$group.labels==unique(env$group.labels)[i])
        group.centroid <- colMeans(secLvlSom$feature.coords[group.member, 1:2])
        
        hull <- chull(secLvlSom$feature.coords[group.member, 1],
                      secLvlSom$feature.coords[group.member, 2])
        
        polygon(secLvlSom$feature.coords[group.member[hull], 1],
                secLvlSom$feature.coords[group.member[hull], 2],
                col=transparent.group.colors[i], lty=1,
                border=env$groupwise.group.colors[i])
      }
      
      legend("topright", as.character(unique(env$group.labels)), cex=0.5,
             text.col=env$groupwise.group.colors, bg="white")
      
      box()
    }
    
    
    ##### Plot SampleSOM with expression portraits ######
    par(mar=c(1,1,1,1))
    xl <- c(0,21)
    yl <- c(0,21)
  
    plot( 0, type="n", axes=FALSE, xlab="", ylab="", xlim=xl, ylim=yl, xaxs="i", yaxs="i",main=paste("Sample SOM on",n), cex.main=1)  
    
    if (length(unique(env$group.labels)) > 1)
    {
      legend("topright", as.character(unique(env$group.labels)), cex=0.5,
             text.col=env$groupwise.group.colors, bg="white")
    }
  
    coord.bins <- as.numeric( cut( c(1:env$preferences$dim.2ndLvlSom), breaks=20 ) )
    for (j in 1:nrow(secLvlSom$node.summary))
    {
      
      which.samples <-
        intersect(which(coord.bins[secLvlSom$feature.coords[,"x"]+1] == coord.bins[secLvlSom$node.summary[j,"x"]+1]),
                  which(coord.bins[secLvlSom$feature.coords[,"y"]+1] == coord.bins[secLvlSom$node.summary[j,"y"]+1]))
  
      if (!is.na(which.samples[1]))
      {
        m <- matrix(env$metadata[, which.samples[1]], env$preferences$dim.1stLvlSom, env$preferences$dim.1stLvlSom)
  
        if (max(m) - min(m) != 0)
        {
          m <- (m - min(m)) / (max(m) - min(m)) * 999
        }
  
        m <- cbind(apply(m, 1, function(x){x}))[nrow(m):1,]
        x <- pixmapIndexed(m , col=env$color.palette.portraits(1000), cellres=10)
  
        addlogo(x,
                coord.bins[secLvlSom$feature.coords[which.samples[1],"x"]+1]+c(-0.45,0.455),
                coord.bins[secLvlSom$feature.coords[which.samples[1],"y"]+1]+c(-0.45,0.45))
  
        which.samples <- which.samples[1:min(9, length(which.samples))]
  
        x.seq <- c(0,0.3,0,-0.3,0,-0.3,0.3,-0.3,0.3)[seq_along(which.samples)]
        y.seq <- c(0,0,0.3,0,-0.3,-0.3,-0.3,0.3,0.3)[seq_along(which.samples)]
  
        points(coord.bins[secLvlSom$feature.coords[which.samples[1],"x"]+1]+x.seq,
               coord.bins[secLvlSom$feature.coords[which.samples[1],"y"]+1]+y.seq,
               pch=16, col=env$group.colors[which.samples], cex=1.2)
  
        points(coord.bins[secLvlSom$feature.coords[which.samples[1],"x"]+1]+x.seq,
               coord.bins[secLvlSom$feature.coords[which.samples[1],"y"]+1]+y.seq,
               pch=1, col="gray20", cex=1.2, lwd=1)
      }
    }
  
    box()
    
  }  
    
  dev.off()
}
