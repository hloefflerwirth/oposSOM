pipeline.sampleSimilarityAnalysisICA <- function(env)
{
  filename <- file.path("Sample Similarity Analysis", "Component Analysis.pdf")
  util.info("Writing:", filename)

  pdf(filename, 21/2.54, 29.7/2.54, useDingbats=FALSE)

  for (i in 1:2 )
  {
    d <- if( i == 1 ) env[[paste("spot.list.",env$preferences$standard.spot.modules,sep="")]]$spotdata else env$metadata
    n <- if( i == 1 ) paste( "module data (", env$preferences$standard.spot.modules, ")") else "metagene data"
    
    try({
      suppressMessages({ ICA.res <- fastICA(t(d), 3)$S })
      z <- ICA.res[,2]

      layout(matrix(c(0, 1, 0)),heights=c(1, 3, 1))
      par(mar=c(1, 1, 1, 1))

      scatterplot3d(ICA.res,
                    cex.symbols=3*(1-((z-min(z))/(max(z)-min(z))))+1,
                    color=env$group.colors,
                    pch=16, xlab="", ylab="", zlab="",
                    main=paste("Independent Component Analysis on",n), mar=c(1,1,1,1))

      par(new=TRUE)
      plot(0,type="n",axes=FALSE,xlab="",ylab="")
      text(0.75, 0.7485, "Component 2", cex=1, srt=38)
      mtext("component 1", 1, cex=0.8, line=-1, at=0.84)
      mtext("component 3", 2, cex=0.8, line=-1, at=-0.3)

      par(new=TRUE)
      plot(0,type="n", axes=FALSE, xlab="", ylab="")

      legend("bottomright", as.character(unique(env$group.labels)), cex=0.5,
             text.col=env$groupwise.group.colors, bg="white")

      
      layout(matrix(c(1,2)))
  
      if(ncol(ICA.res)==3)
      {    
        par(mar=c(0.1,3,1,3))
        
        plot(ICA.res[,1], ICA.res[,3], type="p", pch=16,
             col=env$group.colors, cex=1, axes=FALSE, xlab="", ylab="",
             main=paste("Independent Component Analysis on",n), cex.main=0.8)
        
        mtext("component 3",2,cex=0.8)
        # points(ICA.res[,1], ICA.res[,3], pch=1, col="black", cex=3)
        box()
      }
      
      
      par(mar=c(1,3,0.1,3))
      
      plot(ICA.res[,1], ICA.res[,2], type="p", pch=16,
           col=env$group.colors, cex=1, axes=FALSE, xlab="", ylab="", main="")
      
      mtext("component 1",1,cex=0.8)
      mtext("component 2",2,cex=0.8)
      
      # points(ICA.res[,1], ICA.res[,2], pch=1, col="black", cex=3)
      box()
        
        
      if( ncol(d) < 1000 )
      {
        layout(matrix(c(1,2)))
        
        if(ncol(ICA.res)==3)
        { 
          par(mar=c(0.1,3,1,3))
    
          plot(ICA.res[,1], ICA.res[,3], type="p", pch=16,
               col=env$group.colors, cex=1, axes=FALSE, xlab="", ylab="",
               main=paste("Independent Component Analysis on",n), cex.main=0.8)
    
          mtext("component 3",2,cex=0.8)
          # points(ICA.res[,1], ICA.res[,3], pch=1, col="black", cex=3)
          text(ICA.res[,1], ICA.res[,3], colnames(env$indata), col="gray20", cex=0.6)
          box()
        }
        
        par(mar=c(1,3,0.1,3))
  
        plot(ICA.res[,1], ICA.res[,2], type="p", pch=16,
             col=env$group.colors, cex=1, axes=FALSE, xlab="", ylab="", main="")
  
        mtext("component 1",1,cex=0.8)
        mtext("component 2",2,cex=0.8)
        
        # points(ICA.res[,1], ICA.res[,2], pch=1, col="black", cex=3)
        text(ICA.res[,1], ICA.res[,2], colnames(env$indata), col="gray20", cex=0.6)
        box()
        
        if (length(unique(env$group.labels)) > 1)
        {
          legend("topright", as.character(unique(env$group.labels)), cex=0.4,
                 text.col=env$groupwise.group.colors, bg="white")
        }
      }
      
    }, silent=TRUE)
  }

  dev.off()
}
