pipeline.sampleSimilarityAnalysisCor <- function()
{
  mixColors<-function(c1, c2,alpha=122 )
  {
    do.call( rgb,
             as.list( c( colMeans( sapply(c(2,4,6),function(i)
               c( strtoi(substr(c1,i,i+1),16), strtoi(substr(c2,i,i+1),16) )
               
             ) ), alpha ) /255 )
    )
  }
  
  
  groupwise.group.colors.stability.1 <- apply((col2rgb(groupwise.group.colors) + 0.7 *
           (255 - col2rgb(groupwise.group.colors))) / 255,
        2, function(x) rgb(x[1],x[2],x[3]))
  
  groupwise.group.colors.stability.2 <- apply((col2rgb(groupwise.group.colors) + 0.9 *
          (255 - col2rgb(groupwise.group.colors))) / 255,
       2, function(x) rgb(x[1],x[2],x[3]))  
  
  
  stability.colors <- group.colors
  stability.colors[which(group.silhouette.coef < 0.5)] <- groupwise.group.colors.stability.1[group.labels[names(which(group.silhouette.coef < 0.5))]]
  stability.colors[which(group.silhouette.coef < 0)] <- groupwise.group.colors.stability.2[group.labels[names(which(group.silhouette.coef < 0))]]
    
  pcm.module <- cor( get(paste("spot.list.",preferences$standard.spot.modules,sep=""))$spotdata )
  pcm.metadata <- cor( metadata )
  
  if( any(is.na(pcm.module)) )
  {
    util.warn("Constant sample values in module data. Correlation values of those samples cannot be calculated.")
  }
    
  #### Correlation Spanning Tree ####
  
  filename <- file.path(paste(files.name, "- Results"), "Sample Similarity Analysis", "Correlation Spanning Tree.pdf")
  util.info("Writing:", filename)
  pdf(filename, 21/2.54, 21/2.54, useDingbats=FALSE)
  
  for (i in 1:2 )
  {
    adj.matrix <- if( i == 1 ) pcm.module else pcm.metadata
    adj.matrix[ which(is.na(adj.matrix)) ] <- 0
    n <- if( i == 1 ) paste( "module data (", preferences$standard.spot.modules, ")") else "metagene data"
    
    g <- graph.adjacency(-adj.matrix, weighted=TRUE, mode="undirected")
    stg <- minimum.spanning.tree(g)
    E(stg)$weight <- 1
    layout <- layout_with_kk(stg)
    
    par(mar=c(1,1,1,1))
    E(stg)$color <- apply( get.edgelist( stg ), 1, function(x) mixColors( group.colors[x[1]], group.colors[x[2]] ) )
    plot(stg, layout=layout, vertex.size=5, vertex.label = rep("",ncol(adj.matrix)),
         vertex.color=group.colors, main=paste("Correlation Spanning Tree on",n))
      legend("bottomright", as.character(unique(group.labels)), cex=0.5, text.col=groupwise.group.colors, bg="white")
      box()
 
    E(g)$color <- "darkgrey"
    plot(stg, layout=layout, vertex.size=5, vertex.label = rep("",ncol(adj.matrix)),
             vertex.color=stability.colors, main=paste("Correlation Spanning Tree & silhouette scores on",n))
      legend("bottomright", rep(as.character(unique(group.labels)),3), cex=0.5,
             text.col=c(groupwise.group.colors,groupwise.group.colors.stability.1,groupwise.group.colors.stability.2), bg="white",ncol=3,
             title="0.5<S    0<S<0.5    S<0",title.col="black")
      box()

    if (ncol(adj.matrix) < 1000)
    {
      plot(stg, layout=layout, vertex.size=5, vertex.label = colnames(adj.matrix),
           vertex.label.cex=if (ncol(adj.matrix)<100) 1.2 else 0.6,
           vertex.color=group.colors, main=paste("Correlation Spanning Tree on",n))
        legend("bottomright", as.character(unique(group.labels)), cex=0.5, text.col=groupwise.group.colors, bg="white")
        box()
    }
  }
  dev.off()
  
  #### Correlation Backbone ####
  
  filename <- file.path(paste(files.name, "- Results"), "Sample Similarity Analysis", "Correlation Backbone.pdf")
  util.info("Writing:", filename)
  pdf(filename, 21/2.54, 21/2.54, useDingbats=FALSE)
  
  for (i in 1:2 )
  {
    adj.matrix <- if( i == 1 ) pcm.module else pcm.metadata
    adj.matrix[ which(is.na(adj.matrix)) ] <- 0
    n <- if( i == 1 ) paste( "module data (", preferences$standard.spot.modules, ")") else "metagene data"
    
    diag(adj.matrix) <- 0
    adj.matrix <- apply(adj.matrix, 1, function(x)
    {
      x[order(x,decreasing=TRUE)[-c(1:2)]] <- 0
      return(x)
    })
    adj.matrix[which(adj.matrix < 0.5)] <- 0

    if (max(adj.matrix) > 0)
    {
      g <- graph.adjacency(adj.matrix, weighted=TRUE,  mode="undirected")
      E(g)$weight <- (2 + E(g)$weight)/2
      layout <- layout_with_kk(g)
      
      par(mar=c(1,1,1,1))
      E(g)$color <- apply( get.edgelist( g ), 1, function(x) mixColors( group.colors[x[1]], group.colors[x[2]] ) )
      plot(g, layout=layout, vertex.size=5, vertex.label = rep("",ncol(adj.matrix)),
          vertex.color=group.colors, main=paste("Correlation backbone (2NN-graph) on",n))
        legend("bottomright", as.character(unique(group.labels)), cex=0.5, text.col=groupwise.group.colors, bg="white")
        box()
        
      E(g)$color <- "darkgrey"
      plot(g, layout=layout, vertex.size=5, vertex.label = rep("",ncol(adj.matrix)),
          vertex.color=stability.colors, main=paste("Correlation backbone & silhouette scores on",n))
       legend("bottomright", rep(as.character(unique(group.labels)),3), cex=0.5,
              text.col=c(groupwise.group.colors,groupwise.group.colors.stability.1,groupwise.group.colors.stability.2), bg="white",ncol=3,
              title="0.5<S    0<S<0.5    S<0",title.col="black")
       box()
  
      if (ncol(adj.matrix) < 1000)
      {
        plot(g, layout=layout, vertex.size=5, vertex.label = colnames(adj.matrix),
            vertex.label.cex=if (ncol(adj.matrix)<100) 1.2 else 0.6,
            vertex.color=group.colors, main=paste("Correlation backbone on",n))
         legend("bottomright", as.character(unique(group.labels)), cex=0.5, text.col=groupwise.group.colors, bg="white")
         box()
      }
    }
  }
  dev.off()
  
  #### Correlation Network ####
  
  filename <- file.path(paste(files.name, "- Results"), "Sample Similarity Analysis", "Correlation Network.pdf")
  util.info("Writing:", filename)
  pdf(filename, 21/2.54, 21/2.54, useDingbats=FALSE)
  
  for (i in 1:2 )
  {
    adj.matrix <- if( i == 1 ) pcm.module else pcm.metadata
    adj.matrix[ which(is.na(adj.matrix)) ] <- 0
    n <- if( i == 1 ) paste( "module data (", preferences$standard.spot.modules, ")") else "metagene data"
   
    diag(adj.matrix) <- 0
    adj.matrix[which(adj.matrix < 0.5)] <- 0
    
    if (max(adj.matrix) > 0)
    {
      g <- graph.adjacency(adj.matrix, weighted=TRUE,  mode="undirected")
      E(g)$weight <- (2 + E(g)$weight)/2
      layout <- layout_with_kk( g )

      par(mar=c(1,1,1,1))      
      E(g)$color <- apply( get.edgelist( g ), 1, function(x) mixColors( group.colors[x[1]], group.colors[x[2]], alpha=40 ) )
      plot(g, layout=layout, vertex.size=ifelse(ncol(adj.matrix)<250, 5, 3),
           vertex.label = rep("",ncol(adj.matrix)),
           vertex.color=group.colors, main=paste("Correlation Network on",n))
        legend("bottomright", as.character(unique(group.labels)), cex=0.5, text.col=groupwise.group.colors, bg="white")
        box()
        
      E(g)$color <- "darkgrey"
      plot(g, layout=layout, vertex.size=ifelse(ncol(adj.matrix) < 250, 5, 3),
           vertex.label = rep("",ncol(adj.matrix)),
           vertex.color=stability.colors, main=paste("Correlation Network & silhouette scores on",n))
        legend("bottomright", rep(as.character(unique(group.labels)),3), cex=0.5,
             text.col=c(groupwise.group.colors,groupwise.group.colors.stability.1,groupwise.group.colors.stability.2), bg="white",ncol=3,
             title="0.5<S    0<S<0.5    S<0",title.col="black")
        box()        

      if (ncol(adj.matrix) < 1000)
      {
        plot(g, layout=layout, vertex.size=ifelse(ncol(adj.matrix) < 250, 5, 3),
             vertex.label=colnames(adj.matrix),
             vertex.label.cex=if (ncol(adj.matrix)<100) 1.2 else 0.6,
             vertex.color=group.colors, main=paste("Correlation Network on",n))
          legend("bottomright", as.character(unique(group.labels)), cex=0.5, text.col=groupwise.group.colors, bg="white")
          box()
      }
    }
  }
  dev.off()
  
  #### Pairwise Correlation Maps ####
  
  filename <- file.path(paste(files.name, "- Results"), "Sample Similarity Analysis", "Correlation Maps.pdf")
  util.info("Writing:", filename)
  pdf(filename, 29.7/2.54, 21/2.54, useDingbats=FALSE)
  
  for (i in 1:2 )
  {
    d <- if( i == 1 ) pcm.module else pcm.metadata
    n <- if( i == 1 ) paste( "module data (", preferences$standard.spot.modules, ")") else "metagene data"
    
    d.noNA <- d
    d.noNA[ which(is.na(d.noNA)) ] <- 0
    hcl <- hclust(dist(d.noNA))
    
    par(mar=c(1,1,1,1))
    heatmap(x=d, zlim=c(-1,1), Rowv=NA, Colv=NA, col=color.palette.heatmaps(1000),
                labRow=if(nrow(d)<100) rownames(d) else rep("",nrow(d)),
                labCol=if(ncol(d)<100) colnames(d) else rep("",ncol(d)),
                scale="n", main=paste("Pairwise correlation map on",n),
                margins=c(8,6), ColSideColors=group.colors, RowSideColors=group.colors)
    
    par(new=TRUE, mar=c(5,1,1,2))
    plot(0,type="n", axes=FALSE, xlab="", ylab="")
      legend("bottomright", as.character(unique(group.labels)), cex=0.5, text.col=groupwise.group.colors, bg="white")
      
    par(new=TRUE, mar=c(31.6,55,4.2,2))
    image(matrix(1:100, 1, 100), col=color.palette.heatmaps(1000), axes=FALSE)
       axis(2, c(-1,-0.5,0.5,1,"r"), at=c(0, 0.25, 0.75, 1, 0.5), las=2, tick=FALSE, pos=0, cex.axis=1)

    o <- unlist(sapply(unique(group.labels), function(gr)
    {
      idx <- names(group.labels)[which(group.labels == gr)]

      if (length(idx) > 1)
      {
        hc <- hclust(dist(t(d.noNA[,idx])))
        return(hc$labels[hc$order])
      }
      return(idx)
    }))

    par(mar=c(1,1,1,1))
    heatmap(x=d[o,o], zlim=c(-1,1), Rowv=NA, Colv=NA, col=color.palette.heatmaps(1000),
                  labRow=if(nrow(d)<100) rownames(d[o]) else rep("",nrow(d)),
                  labCol=if(ncol(d)<100) colnames(d[o]) else rep("",ncol(d)),
                  scale="n", main=paste("Pairwise correlation map on",n),
                  margins=c(8,6), ColSideColors=group.colors[o], RowSideColors=group.colors[o])

    par(new=TRUE, mar=c(5,1,1,2))
    plot(0,type="n", axes=FALSE, xlab="", ylab="")
    legend("bottomright", as.character(unique(group.labels)), cex=0.5, text.col=groupwise.group.colors, bg="white")

    par(new=TRUE, mar = c(31.6, 55, 4.2, 2))
    image(matrix(1:100, 1, 100), col=color.palette.heatmaps(1000), axes=FALSE)
       axis(2, c(-1,-0.5,0.5,1,"r"), at=c(0, 0.25, 0.75, 1, 0.5), las=2, tick=FALSE, pos=0, cex.axis=1)

       
    par(mar=c(1,1,1,1))
    heatmap(x=d, zlim=c(-1,1), Rowv=as.dendrogram(hcl), Colv=as.dendrogram(hcl),
                 labRow=if(nrow(d)<100) rownames(d) else rep("",nrow(d)),
                 labCol=if(ncol(d)<100) colnames(d) else rep("",ncol(d)),
                 col=color.palette.heatmaps(1000),
                 scale="n", main=paste("Pairwise correlation map on",n),
                 margins=c(8,6), ColSideColors=group.colors, RowSideColors=group.colors)

    par(new=TRUE, mar=c(5,1,1,2))
    plot(0,type="n", axes=FALSE, xlab="", ylab="")
    legend("bottomright", as.character(unique(group.labels)), cex=0.5,
           text.col=groupwise.group.colors, bg="white")

    par(new=TRUE, mar = c(25, 55, 10.8, 2))
    image(matrix(1:100, 1, 100), col=color.palette.heatmaps(1000), axes=FALSE)
      axis(2, c(-1,-0.5,0.5,1,"r"), at=c(0, 0.25, 0.75, 1, 0.5), las=2, tick=FALSE, pos=0, cex.axis=1)

    if (ncol(d) < 1000 && i == 1)
    {
      def.par <- par(no.readonly = TRUE)
      
      environment(pipeline.moduleCorrelationMap) <- environment()
      pipeline.moduleCorrelationMap(d.noNA, hcl)
      
      par(def.par)
    }
    
  }
  dev.off()

}
