pipeline.3rdLvlNetworks <- function()
{
  plot.set.list.networks <- function(set.list, main)
  {
    N.spots <- length(set.list$spots)

    if (N.spots > 1)
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

      environment(pipeline.3rdLvlNetworksWto) <- environment()
      pipeline.3rdLvlNetworksWto(set.list, sample.spot.matrix)

      # plot single spot MST
      mean.spot.expression <- matrix(NA, ncol(indata), 0)

      for (i in 1:length(set.list$spots))
      {
        spot.metagenes <- set.list$spots[[i]]$metagenes

        mean.spot.expression <- cbind(mean.spot.expression,
                                      if (length(spot.metagenes) > 1)
                                        colMeans(metadata[spot.metagenes,])
                                      else
                                        metadata[spot.metagenes,])
      }

      adj.matrix <- cor(mean.spot.expression)
      diag(adj.matrix) <- 0
      adj.matrix[which(adj.matrix < 0.5)] <- 0

      if (max(adj.matrix) > 0)
      {
        g <- graph.adjacency(adj.matrix, weighted=T,  mode="undirected")

        layout <- layout.fruchterman.reingold(g)
        layout <- layout.norm(layout, xmin=-1, xmax=1, ymin=-1, ymax=1)

        par(mar=c(1,1,1,1), mfrow=c(1,1))

        plot(g, layout=layout, vertex.size=0, vertex.label="",
             xlim=c(-1.1,1.1), ylim=c(-1.1,1.1))

        text(-1.65, 1.1, "Single spot network", adj=0, cex=2)

        w <- 0.08

        for (i in 1:N.spots)
        {
          samples <- which(sample.spot.matrix[i,] == 1)
          tbl <- table(as.vector(group.labels[samples]))

          for (ii in 1:length(tbl))
          {
            tbl[ii] <- tbl[ii] / table(as.vector(group.labels))[names(tbl[ii])]
          }

          tbl <- round(100 * tbl)
          t2 <- table(as.vector(group.labels[samples]))

          mask <- matrix(set.list$spots[[i]]$mask, preferences$dim.1stLvlSom, preferences$dim.1stLvlSom)
          mask <- cbind(apply(mask, 1, function(x){x}))[nrow(mask):1,]
          mask[is.na(mask)] <- 0

          rect(layout[i,1] - w,
               layout[i,2] - w,
               layout[i,1] + w,
               layout[i,2] + w,
               col="white")

          x <- pixmapIndexed(mask , col = c("white", "green"))
          addlogo(x, layout[i,1] + c(-w,w), layout[i,2] + c(-w,w))

          rect(layout[i,1] - w,
               layout[i,2] - w,
               layout[i,1] + w,
               layout[i,2] + w)

          text(layout[i,1]-w+0.01, layout[i,2]+w+0.025,
               paste("Total :", length(samples)), adj=0, cex=0.6)

          for (ii in 1:length(tbl))
          {
            text(layout[i,1]-w+0.01, layout[i,2]+w+0.01  - ii*0.036,
                 paste(names(tbl)[ii], ":", t2[ii], " - ", tbl[ii], "%"),
                 adj=0, cex=0.6, col=group.colors[which(group.labels == names(tbl)[ii])[1]])
          }
        }
      }

      adj.matrix <- cor(mean.spot.expression) * -1

      g <- graph.adjacency(adj.matrix, weighted=T, mode="undirected")
      stg <- minimum.spanning.tree(g)

      layout <- layout.kamada.kawai(stg)
      layout <- layout.norm(layout, xmin=-1, xmax=1, ymin=-1, ymax=1)

      par(mar=c(1,1,1,1), mfrow=c(1,1))
      plot(stg, layout=layout, vertex.size=0, vertex.label="", xlim=c(-1.1,1.1), ylim=c(-1.1,1.1))
      text(-1.65, 1.1, "Correlation Spanning Tree", adj=0, cex=2)
      text(-1.65, 1, "Single Spots", adj=0, cex=2)

      w <- 0.08

      for (i in 1:N.spots)
      {
        samples <- which(sample.spot.matrix[i,] == 1)
        tbl <- table(as.vector(group.labels[samples]))

        for (ii in 1:length(tbl))
        {
          tbl[ii] <- tbl[ii] / table(as.vector(group.labels))[names(tbl[ii])]
        }

        tbl <- round(100 * tbl)
        t2 <- table(as.vector(group.labels[samples]))

        mask <- matrix(set.list$spots[[i]]$mask, preferences$dim.1stLvlSom, preferences$dim.1stLvlSom)
        mask <- cbind(apply(mask, 1, function(x){x}))[nrow(mask):1,]
        mask[is.na(mask)] <- 0

        rect(layout[i,1] - w,
             layout[i,2] - w,
             layout[i,1] + w,
             layout[i,2] + w,
             col="white")

        x <- pixmapIndexed(mask , col = c("white", "green"))
        addlogo(x, layout[i,1] + c(-w,w), layout[i,2] + c(-w,w))

        rect(layout[i,1] - w,
             layout[i,2] - w,
             layout[i,1] + w,
             layout[i,2] + w)

        text(layout[i,1]-w+0.01, layout[i,2]+w+0.025,
             paste("Total :", length(samples)), adj=0, cex=0.6)

        for (ii in 1:length(tbl))
        {
          text(layout[i,1]-w+0.01, layout[i,2]+w+0.01  - ii*0.036,
               paste(names(tbl)[ii], ":", t2[ii], " - ", tbl[ii], "%"),
               adj=0, cex=0.6, col=group.colors[which(group.labels == names(tbl)[ii])[1]])
        }
      }

      if (length(unique(group.labels)) > 1)
      {
        plot(stg, layout=layout, vertex.size=0, vertex.label="",
             xlim=c(-1.1,1.1), ylim=c(-1.1,1.1))

        text(-1.65, 1.1, "Correlation Spanning Tree", adj=0, cex=2)
        text(-1.65, 1, "absolute numbers", adj=0)

        par(new=T)

        tbl <- matrix(0, N.spots, length(unique(group.labels)))
        colnames(tbl) <- unique(group.labels)

        for (i in 1:N.spots)
        {
          samples <- which(sample.spot.matrix[i,] == 1)

          tbl[i, names(table(as.vector(group.labels[samples])))] <-
            table(as.vector(group.labels[samples]))
        }

        tbl <- tbl / max(tbl)

        stars(tbl, locations=layout, xlim=c(-1.1,1.1), ylim=c(-1.1,1.1),
              len=0.1, draw.segments=T, col.segments=groupwise.group.colors,
              key.loc=c(1.3,-1.1), scale=F)

        plot(stg, layout=layout, vertex.size=0, vertex.label="",
             xlim=c(-1.1,1.1), ylim=c(-1.1,1.1))

        text(-1.65, 1.1, "Correlation Spanning Tree", adj=0, cex=2)
        text(-1.65, 1, "relative to category size", adj=0)

        par(new=T)

        tbl <- matrix(0, N.spots, length(unique(group.labels)))
        colnames(tbl) <- unique(group.labels)

        for (i in 1:N.spots)
        {
          samples <- which(sample.spot.matrix[i,] == 1)
          tt <- table(as.vector(group.labels[samples]))

          for (ii in 1:length(tt))
          {
            tt[ii] <- tt[ii] / table(as.vector(group.labels))[names(tt[ii])]
          }

          tbl[i, names(tt)] <- tt
        }

        stars(tbl, locations=layout, xlim=c(-1.1,1.1), ylim=c(-1.1,1.1),
              len=0.1, draw.segments=T, col.segments=groupwise.group.colors,
              key.loc=c(1.3,-1.1), scale=F)
      }
    }
  }

  dirname <- file.path(paste(files.name, "- Results"), "3rd lvl Spot Analysis")

  filename <-
    if (spot.list.overexpression$filtered)
    {
      file.path(dirname, "Overexpression Networks filtered.pdf")
    } else
    {
      file.path(dirname, "Overexpression Networks.pdf")
    }

  util.info("Writing:", filename)
  pdf(filename, 29.7/2.54, 21/2.54)
  plot.set.list.networks(set.list=spot.list.overexpression, main="Sample-Overexpression")
  dev.off()

  filename <-
    if (spot.list.underexpression$filtered)
    {
      file.path(dirname, "Underexpression Networks filtered.pdf")
    } else
    {
      file.path(dirname, "Underexpression Networks.pdf")
    }

  util.info("Writing:", filename)
  pdf(filename, 29.7/2.54, 21/2.54)
  plot.set.list.networks(set.list=spot.list.underexpression, main="Sample-Underexpression")
  dev.off()

  filename <- file.path(dirname, "K-Means Cluster Networks.pdf")
  util.info("Writing:", filename)
  pdf(filename, 29.7/2.54, 21/2.54)
  plot.set.list.networks(set.list=spot.list.kmeans, main="K-Means Cluster")
  dev.off()

  if (length(unique(group.labels)) > 1)
  {
    filename <- file.path(dirname, "Group Overexpression Networks.pdf")
    util.info("Writing:", filename)
    pdf(filename, 29.7/2.54, 21/2.54)
    plot.set.list.networks(set.list=spot.list.group.overexpression, main="Group Overexpression")
    dev.off()
  }
}
