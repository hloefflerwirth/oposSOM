pipeline.sampleSimilarityAnalysisED <- function()
{

  #### t-SNE ####
  
  filename <- file.path(paste(files.name, "- Results"), "Sample Similarity Analysis", "t-SNE.pdf")
  util.info("Writing:", filename)
  pdf(filename, 21/2.54, 21/2.54, useDingbats=FALSE)

  for (i in 1:2 )
  {
    d <- if( i == 1 ) get(paste("spot.list.",preferences$standard.spot.modules,sep=""))$spotdata else metadata
    n <- if( i == 1 ) paste( "module data (", preferences$standard.spot.modules, ")") else "metagene data"
    
    suppressMessages({  tsne.res <- tsne( t(d), perplexity=max(30,ncol(indata)/20), max_iter=3000 )  })
    
    par(mar=c(2,2,2,2))
    plot(tsne.res, col=group.colors, pch=16)
      title(paste("t-SNE on",n))  
      legend("bottomright", as.character(unique(group.labels)), cex=0.5, text.col=groupwise.group.colors, bg="white")
  }
  dev.off()
  
  #### Neighbor Joining ####
  
  filename <- file.path(paste(files.name, "- Results"), "Sample Similarity Analysis", "Neighbor Joining.pdf")
  util.info("Writing:", filename)
  pdf(filename, 21/2.54, 21/2.54, useDingbats=FALSE)
  
  for (i in 1:2 )
  {
    d <- if( i == 1 ) get(paste("spot.list.",preferences$standard.spot.modules,sep=""))$spotdata else metadata
    n <- if( i == 1 ) paste( "module data (", preferences$standard.spot.modules, ")") else "metagene data"

    phylo.tree <- nj(dist(t(d)))
    phylo.tree$tip.label <- colnames(d)
    
    par(mar=c(1,1,1,1))
    plot.phylo(phylo.tree, "unrooted", cex=0.5, tip.color=group.colors)
      title(paste("Phylogenetic clustering on",n))
      box()
    
    par(new=TRUE)
    plot(0,type="n", axes=FALSE, xlab="", ylab="")    
    legend("bottomright", as.character(unique(group.labels)), cex=0.5, text.col=groupwise.group.colors, bg="white")
  }
  dev.off()

  #### Heatmaps ####

  filename <- file.path(paste(files.name, "- Results"), "Sample Similarity Analysis", "Hierarchical Clustering.pdf")
  util.info("Writing:", filename)
  pdf(filename, 29.7/2.54, 21/2.54, useDingbats=FALSE)
  
  for (i in 1:2 )
  {
    d <- if( i == 1 ) get(paste("spot.list.",preferences$standard.spot.modules,sep=""))$spotdata else metadata
    n <- if( i == 1 ) paste( "module data (", preferences$standard.spot.modules, ")") else "metagene data"
    
    par(mar=c(1,1,1,1))
    heatmap(x=d, col=color.palette.heatmaps(1000), main=paste("Heatmap on",n), scale="n",
                 labCol=if(ncol(d)<100) colnames(d) else rep("",ncol(d)),
                 labRow=if(i==1) NULL else NA,
                 margins=c(8, 6), zlim=max(abs(range(d)))*c(-1,1), ColSideColors=group.colors, Colv=NA)
    
    par(new=TRUE, mar=c(5,1,1,2))
    plot(0,type="n", axes=FALSE, xlab="", ylab="")
      legend("bottomright", as.character(unique(group.labels)), cex=0.5, text.col=groupwise.group.colors, bg="white")
      
    par(new=TRUE, mar = c(31.6, 55, 4.2, 2))
    image(matrix(1:100, 1, 100), col = color.palette.heatmaps(1000), axes=FALSE)
      axis(2, round(max(abs(range(d)))*c(-1,1),1), at=c(0, 1), las=2, tick=FALSE, pos=0, cex.axis=1)
    
    par(mar=c(1,1,1,1))
    heatmap(x=d, col=color.palette.heatmaps(1000), main=paste("Clustering heatmap on",n), scale="n",
                 labCol=if(ncol(d)<100) colnames(d) else rep("",ncol(d)),
                 labRow=if(i==1) NULL else NA,
                 margins=c(8, 5), zlim=max(abs(range(d)))*c(-1,1), ColSideColors=group.colors)
    
    par(new=TRUE, mar=c(5,1,1,2))
    plot(0,type="n", axes=FALSE, xlab="", ylab="")
      legend("bottomright", as.character(unique(group.labels)), cex=0.5, text.col=groupwise.group.colors, bg="white")
      
    par(new=TRUE, mar = c(25, 55, 10.8, 2))
    image(matrix(1:100, 1, 100), col = color.palette.heatmaps(1000), axes=FALSE)
      axis(2, round(max(abs(range(d)))*c(-1,1),1), at=c(0, 1), las=2, tick=FALSE, pos=0, cex.axis=1)
  }
  dev.off()

}
