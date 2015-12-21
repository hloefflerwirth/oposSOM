pipeline.3rdLvlNetworks <- function()
{
  plot.set.list.networks <- function(set.list, main)
  {
    N.spots <- length(set.list$spots)

    if (N.spots > 2)
    {
      thresh.global <- sd(as.vector(set.list$spotdata))

      sample.spot.matrix <- matrix(0, length(set.list$spots), ncol(metadata))
      rownames(sample.spot.matrix) <- names(set.list$spots)

      for (s in names(set.list$spots))
      {
        if (main %in% c("Sample-Underexpression"))
        {
          sample.spot.matrix[s,] <-
            as.numeric(set.list$spotdata[s,] < mean(set.list$spotdata[s,])-thresh.global)
        } else
        {
          sample.spot.matrix[s,] <-
            as.numeric(set.list$spotdata[s,] > mean(set.list$spotdata[s,])+thresh.global)
        }
      }

    
      ### calculate wTO network ###
      adj_matrix <- cor(t(set.list$spotdata),t(set.list$spotdata))
      diag(adj_matrix) <- 0
      
      omega <- matrix(NA,
                      nrow=nrow(set.list$spotdata),
                      ncol=nrow(set.list$spotdata),
                      dimnames=list(names(set.list$spots),names(set.list$spots)))
      
      b <- adj_matrix %*% adj_matrix
      k_sum <- apply(adj_matrix, 2, function(x) { sum(abs(x)) })
      
      for (i in 1: nrow(set.list$spotdata))
      {
        for (j in 1: nrow(set.list$spotdata))
        {
          omega[i,j] <-
            (b[i,j] + adj_matrix[i,j]) /
            (min(k_sum[i],k_sum[j]) + 1 - abs(adj_matrix[i,j]))
        }
      }
      
      diag(omega) <- 0
      omega[which(abs(omega) < 0.25)] <- 0
      
      g <- graph.adjacency(omega, weighted=TRUE, mode="undirected")

      V(g)$label <- names(set.list$spots)
      n.spot <- apply(sample.spot.matrix, 1, sum)
      
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
      par(mar=c(4.9,13.5,4.9,13.5))
      image(matrix(set.list$overview.mask, preferences$dim.1stLvlSom), col="gray90", axes=FALSE, main="WTO network mapped to SOM space", cex.main=1.6)
      par(new=TRUE, mar=c(1,1,1,1))
      
      plot(g, layout=t(sapply(set.list$spots, function(x) { x$position })),
           vertex.label.color="black" ,vertex.label.cex=1, vertex.color="grey",
           xlim=c(-1.2,1.2), ylim=c(-1.2,1.2))
      
      rect(-1.1,-1.1,1.1,1.1)
      
      legend("bottomright", c("negative correlation",expression(paste("  ",omega," < -0.25")),expression(paste("  ",omega," < -0.50")),expression(paste("  ",omega," < -0.75")),
                              "positive correlation",expression(paste("  ",omega," > 0.25")),expression(paste("  ",omega," > 0.50")),expression(paste("  ",omega," > 0.75"))),
             lty = c(1,1), lwd = 4,col=c("white","rosybrown1","chocolate","brown4","white","lightblue","cornflowerblue","dodgerblue4"))
    }
  }

  dirname <- file.path(paste(files.name, "- Results"), "3rd lvl Spot Analysis")

  filename <-
    if (spot.list.overexpression$filtered)
    {
      file.path(dirname, "wTO Networks - Overexpression Spots filtered.pdf")
    } else
    {
      file.path(dirname, "wTO Networks - Overexpression Spots.pdf")
    }

  util.info("Writing:", filename)
  pdf(filename, 29.7/2.54, 21/2.54)
  plot.set.list.networks(set.list=spot.list.overexpression, main="Overexpression Spots")
  dev.off()

  filename <-
    if (spot.list.underexpression$filtered)
    {
      file.path(dirname, "wTO Networks - Underexpression Spots filtered.pdf")
    } else
    {
      file.path(dirname, "wTO Networks - Underexpression Spots.pdf")
    }

  util.info("Writing:", filename)
  pdf(filename, 29.7/2.54, 21/2.54)
  plot.set.list.networks(set.list=spot.list.underexpression, main="Underexpression Spots")
  dev.off()

  filename <- file.path(dirname, "wTO Networks - K-Means Clusters.pdf")
  util.info("Writing:", filename)
  pdf(filename, 29.7/2.54, 21/2.54)
  plot.set.list.networks(set.list=spot.list.kmeans, main="K-Means Clusters")
  dev.off()

  if (length(unique(group.labels)) > 1)
  {
    filename <- file.path(dirname, "wTO Networks - Group Overexpression Spots.pdf")
    util.info("Writing:", filename)
    pdf(filename, 29.7/2.54, 21/2.54)
    plot.set.list.networks(set.list=spot.list.group.overexpression, main="Group Overexpression Spots")
    dev.off()
  }
  
  filename <- file.path(dirname, "wTO Networks - D-Clusters.pdf")
  util.info("Writing:", filename)
  pdf(filename, 29.7/2.54, 21/2.54)
  plot.set.list.networks(set.list=spot.list.dmap, main="D-Clusters")
  dev.off()
}
