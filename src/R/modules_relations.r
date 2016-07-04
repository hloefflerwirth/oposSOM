modules.relations <- function(spot.list, main, path)
{
  sd.theshold = sd(spot.list$spotdata)
  
  pdf(path, 29.7/2.54, 21/2.54)

  #### wTO network ####
  
  if ( length(spot.list$spots) > 2 )
  {
    # calculate wTO #
    adj_matrix <- cor(t(spot.list$spotdata),t(spot.list$spotdata))
    diag(adj_matrix) <- 0
    
    omega <- matrix(NA,
                    nrow=nrow(spot.list$spotdata),
                    ncol=nrow(spot.list$spotdata),
                    dimnames=list(names(spot.list$spots),names(spot.list$spots)))
    
    b <- adj_matrix %*% adj_matrix
    k_sum <- apply(adj_matrix, 2, function(x) { sum(abs(x)) })
    
    for (i in 1: nrow(spot.list$spotdata))
    {
      for (j in 1: nrow(spot.list$spotdata))
      {
        omega[i,j] <-
          (b[i,j] + adj_matrix[i,j]) /
          (min(k_sum[i],k_sum[j]) + 1 - abs(adj_matrix[i,j]))
      }
    }
    
    diag(omega) <- 0
    omega[which(abs(omega) < 0.25)] <- 0
    
    g <- graph.adjacency(omega, weighted=TRUE, mode="undirected")

    V(g)$label <- names(spot.list$spots)
    n.spot <- apply(spot.list$spotdata, 1, function(x) sum(x>sd.theshold) )
    
    if (max(n.spot) == min(n.spot))
    {
      V(g)$size <- 10
    } else
    {
      V(g)$size <- 8 * (n.spot-min(n.spot)) / (max(n.spot)-min(n.spot)) + 6
    }
    
    if (length(E(g)) > 0)
    {
      E(g)$color <- c("brown4","chocolate","rosybrown1","white","white","lightblue","cornflowerblue","dodgerblue4")[ cut( E(g)$weight, breaks=c(-Inf,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,Inf), labels=c(1:8) ) ]
      E(g)$width <- abs(E(g)$weight)
      
      if (max(E(g)$width) == min(E(g)$width))
      {
        E(g)$width <- 3
      } else
      {
        E(g)$width <- 5 * (E(g)$width-min(E(g)$width)) /
          (max(E(g)$width)-min(E(g)$width)) + 0.5
      }
    }
    
    # plot wTO network
    
    layout(matrix(c(1, 2), 1, 2), c(2, 1), 1)
    par(mar=c(5, 4, 4, 1.8))
    
    image(matrix(spot.list$overview.mask, preferences$dim.1stLvlSom), col="gray90", axes=FALSE, main="WTO network mapped to SOM space", cex.main=1.6)
      box()
    
    par(new=TRUE, mar=c(3.2,2.2,2,0))
    plot(g, layout=t(sapply(spot.list$spots, function(x) { x$position })),
        vertex.label.color="black" ,vertex.label.cex=1, vertex.color="grey",
        xlim=c(-1.0,1.0), ylim=c(-1.0,1.0))
  
    # plot(g, layout=matrix(c(1,1,20,20,1,20,20),ncol=2),
    #      vertex.label.color="black" ,vertex.label.cex=1, vertex.color="grey",
    #      xlim=c(-1.0,1.0), ylim=c(-1.0,1.0))
    
    par(mar=c(5, 1, 4, 2))
    plot(0, type="n", axes=FALSE, xlab="", ylab="", xlim=c(0,1), ylim=c(0,1), xaxs="i", yaxs="i")
    legend("center", c("negative correlation",expression(paste("  ",omega," < -0.25")),expression(paste("  ",omega," < -0.50")),expression(paste("  ",omega," < -0.75")),
                          "positive correlation",expression(paste("  ",omega," > 0.25")),expression(paste("  ",omega," > 0.50")),expression(paste("  ",omega," > 0.75"))),
                          lty = c(1,1), lwd = 4,col=c("white","rosybrown1","chocolate","brown4","white","lightblue","cornflowerblue","dodgerblue4"), bty="n")
    box()
  }
 
    
  #### group association #### 
  
  if (length(unique(group.labels)) > 1)
  {
    spot.group.assoc <- apply(spot.list$spotdata, 1, function(x) tapply(x>sd.theshold,group.labels,sum)[unique(group.labels)] / table(group.labels) )
    spot.goup.mean.number <- tapply( apply(spot.list$spotdata>sd.theshold, 2, sum ), group.labels, mean)[unique(group.labels)]
    
    # barplot
        
    layout(matrix(c(1, 2), 1, 2), c(2, 1), 1)
    par(mar=c(5, 4, 4, 1 ))
    
    barplot(spot.group.assoc[,ncol(spot.group.assoc):1],
            main="Association of the groups to the spots", cex.main=1.5, cex.axis=2, cex.names=1,
            col=groupwise.group.colors, horiz=TRUE, las=1)
    
    par(mar=c(5, 1, 4, 2))
    
    plot(0, main="", cex.main=1, type="n", axes=FALSE, xlab="",
         ylab="", xlim=c(0,1), ylim=c(0,1), xaxs="i", yaxs="i")
    
      legend("center", legend=paste(unique(group.labels), "(", round(spot.goup.mean.number,1), ")"),
           col=groupwise.group.colors, pch=15, pt.cex=1.6, title="< #spots >", bty="n")
      box()
    
    # mapping to SOM
      
    layout(matrix(c(1, 2), 1, 2), c(2, 1), 1)
    par(mar=c(5, 4, 4, 1))
    
    image(matrix(spot.list$overview.mask, preferences$dim.1stLvlSom), col="gray90", axes=FALSE,
          main="Association of the groups to the spots", cex.main=1.5)
    par(new=TRUE)
    plot(0, type="n", xlab="", ylab="", axes=FALSE, xlim=c(0,preferences$dim.1stLvlSom),
         ylim=c(0,preferences$dim.1stLvlSom), xaxs="i", yaxs="i")
    
    box()
    
    for (i in seq_along(spot.list$spots))
    {
      stars(t(spot.group.assoc[,i]), locations=spot.list$spots[[i]]$position,
            len=2, draw.segments=TRUE,  scale=FALSE, add=TRUE, col.segments=groupwise.group.colors)
      
      par(fg="gray", lty=2)
      
      stars(t(c(0.5)), locations=spot.list$spots[[i]]$position, len=2,
            draw.segments=TRUE,  scale=FALSE, add=TRUE, col.segments=NA)
      stars(t(c(1)), locations=spot.list$spots[[i]]$position, len=2,
            draw.segments=TRUE,  scale=FALSE, add=TRUE, col.segments=NA)
      
      par(fg="black", lty="solid")
      
      text(spot.list$spots[[i]]$position[1], spot.list$spots[[i]]$position[2]+preferences$dim.1stLvlSom*0.05,
           names(spot.list$spots)[i], col="gray50")
    }
    
    par(mar=c(5, 1, 4, 2))
    plot(0, type="n", axes=FALSE, xlab="", ylab="", xlim=c(0,1), ylim=c(0,1), xaxs="i", yaxs="i")
      legend("center", unique(group.labels), col = groupwise.group.colors, pch=15, pt.cex=1.6, bty="n")
      box()
  }  
    
    
    
  dev.off()
}
