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

      pattern.spot.matrix <- matrix(NA, N.spots, 0)
      pattern.sample.order <- list()

      sample.spot.matrix2 <- sample.spot.matrix
      colnmes <- 1:ncol(sample.spot.matrix2)

      p <- 1
      while (ncol(as.matrix(sample.spot.matrix2)) > 0)
      {
        if (ncol(as.matrix(sample.spot.matrix2)) > 1)
        {
          identic.pattern <- c(colnmes[1])

          for (j2 in 2:ncol(sample.spot.matrix2))
          {
            if (all(sample.spot.matrix2[,1] == sample.spot.matrix2[,j2]))
            {
              identic.pattern <- c(identic.pattern, colnmes[j2])
            }
          }

          pattern.spot.matrix <- cbind(pattern.spot.matrix, sample.spot.matrix2[,1])
          pattern.sample.order[[p]] <- identic.pattern

          del <- which(colnmes %in% identic.pattern)
          sample.spot.matrix2 <- sample.spot.matrix2[, -del]
          colnmes <- colnmes[-del]

          p <- p+1
        } else
        {
          pattern.spot.matrix <- cbind(pattern.spot.matrix, sample.spot.matrix2)
          pattern.sample.order[[p]] <- colnmes[1]

          sample.spot.matrix2 <- matrix(NA, 0, 0)
        }
      }

      pattern.adj.matrix <- matrix(NA,
                                   ncol(pattern.spot.matrix),
                                   ncol(pattern.spot.matrix))

      for (j in 1:(ncol(pattern.spot.matrix)-1))
      {
        for (j2 in (j+1):ncol(pattern.spot.matrix))
        {
          common.pattern <- which(pattern.spot.matrix[,j] & pattern.spot.matrix[,j2])
          diff.pattern <- which(pattern.spot.matrix[,j] != pattern.spot.matrix[,j2])
          sum.dist <- rep(NA, nrow(pattern.spot.matrix))

          for (dp in diff.pattern)
          {
            sum.dist[dp] <- 0

            for (cp in common.pattern)
            {
              ed <- sqrt(sum((set.list$spots[[cp]]$position - set.list$spots[[dp]]$position) ^ 2))
              sum.dist[dp] <- sum.dist[dp] + ed
            }

            sum.dist[dp] <- sum.dist[dp] / length(common.pattern)
          }

          sum.dist <- sum(sum.dist, na.rm=T)

          pattern.adj.matrix[j, j2] <- sum.dist
          pattern.adj.matrix[j2, j] <- sum.dist
        }
      }

      pattern.adj.matrix[is.na(pattern.adj.matrix)] <- 2 * max(pattern.adj.matrix, na.rm=T)

      environment(pipeline.3rdLvlNetworksWto) <- environment()
      pipeline.3rdLvlNetworksWto(set.list)

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

          mask <- matrix(set.list$spots[[i]]$mask, preferences$dim.som1, preferences$dim.som1)
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

        mask <- matrix(set.list$spots[[i]]$mask, preferences$dim.som1, preferences$dim.som1)
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
              len=0.1, draw.segments=T, col.segments=unique.group.colors,
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
              len=0.1, draw.segments=T, col.segments=unique.group.colors,
              key.loc=c(1.3,-1.1), scale=F)
      }

      # plot spot pattern MST
      if (ncol(pattern.adj.matrix) > 2)
      {
        adj.matrix <- pattern.adj.matrix
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

          text(-1.65, 1.1, "Spot pattern network", adj=0, cex=2)

          for (j in 1:ncol(pattern.spot.matrix))
          {
            samples <- pattern.sample.order[[j]]
            spots <- which(pattern.spot.matrix[,j] == 1)
            mask <- matrix(0, preferences$dim.som1, preferences$dim.som1)

            for (s in spots)
            {
              msk <- set.list$spots[[s]]$mask
              mask[which(!is.na(msk))] <- 1
            }

            mask <- cbind(apply(mask, 1, function(x){x}))[nrow(mask):1,]

            rect(layout[j,1] - w,
                 layout[j,2] - w,
                 layout[j,1] + w,
                 layout[j,2] + w,
                 col="white")

            x <- pixmapIndexed(mask , col = c("white", "green"))
            addlogo(x, layout[j,1] + c(-w,w), layout[j,2] + c(-w,w))

            rect(layout[j,1] - w,
                 layout[j,2] - w,
                 layout[j,1] + w,
                 layout[j,2] + w)

            gr <- as.vector(group.labels[samples])
            tbl <- table(gr)

            for (i in 1:length(tbl))
            {
              tbl[i] <- tbl[i] / table(as.vector(group.labels))[names(tbl[i])]
            }

            tbl <- round(100 * tbl)
            t2 <- table(gr)

            text(layout[j,1]-w+0.01, layout[j,2]+w+0.025,
                 paste("Total :", length(samples)), adj=0, cex=0.6)

            for (i in 1:length(tbl))
            {
              text(layout[j,1]-w+0.01, layout[j,2]+w+0.01  - i*0.036,
                   paste(names(tbl)[i], ":", t2[i], " - ", tbl[i], "%"),
                   adj=0, cex=0.6, col=group.colors[which(group.labels == names(tbl)[i])[1]])
            }
          }
        }

        g <- graph.adjacency(pattern.adj.matrix, weighted=T,  mode="undirected")
        stg <- minimum.spanning.tree(g)

        if (length(E(stg)) > 0)
        {
          layout <- layout.fruchterman.reingold(stg)
          layout <- layout.norm(layout, xmin=-1, xmax=1, ymin=-1, ymax=1)

          par(mar=c(1,1,1,1), mfrow=c(1,1))
          plot(stg, layout=layout, vertex.size=0, vertex.label="", xlim=c(-1.1,1.1), ylim=c(-1.1,1.1))
          text(-1.65, 1.1, "Distance Spanning Tree", adj=0, cex=2)
          text(-1.65, 1, "Spot Pattern", adj=0, cex=2)

          for (j in 1:ncol(pattern.spot.matrix))
          {
            samples <- pattern.sample.order[[j]]
            spots <- which(pattern.spot.matrix[,j] == 1)
            mask <- matrix(0, preferences$dim.som1, preferences$dim.som1)

            for (s in spots)
            {
              msk <- set.list$spots[[s]]$mask
              mask[which(!is.na(msk))] <- 1
            }

            mask <- cbind(apply(mask, 1, function(x){x}))[nrow(mask):1,]
            rect(layout[j,1] - w, layout[j,2] - w, layout[j,1] + w, layout[j,2] + w, col="white")
            x <- pixmapIndexed(mask , col = c("white", "green"))
            addlogo(x, layout[j,1] + c(-w,w), layout[j,2] + c(-w,w))

            rect(layout[j,1] - w,
                 layout[j,2] - w,
                 layout[j,1] + w,
                 layout[j,2] + w)

            gr <- as.vector(group.labels[samples])
            tbl <- table(gr)

            for (i in 1:length(tbl))
            {
              tbl[i] <- tbl[i] / table(as.vector(group.labels))[names(tbl[i])]
            }

            tbl <- round(100 * tbl)
            t2 <- table(gr)

            text(layout[j,1]-w+0.01, layout[j,2]+w+0.025,
                 paste("Total :", length(samples)), adj=0, cex=0.6)

            for (i in 1:length(tbl))
            {
              text(layout[j,1]-w+0.01, layout[j,2]+w+0.01  - i*0.036,
                   paste(names(tbl)[i], ":", t2[i], " - ", tbl[i], "%"),
                   adj=0, cex=0.6, col=group.colors[which(group.labels == names(tbl)[i])[1]])
            }
          }

          if (length(unique(group.labels)) > 1)
          {
            plot(stg, layout=layout, vertex.size=0, vertex.label="",
                 xlim=c(-1.1,1.1), ylim=c(-1.1,1.1))

            text(-1.65, 1.1, "Distance Spanning Tree", adj=0, cex=2)
            text(-1.65, 1, "absolute numbers", adj=0)

            par(new=T)
            tbl <- matrix(0, ncol(pattern.adj.matrix), length(unique(group.labels)))
            colnames(tbl) = unique(group.labels)

            for (i in 1:ncol(pattern.adj.matrix))
            {
              samples <- pattern.sample.order[[i]]

              tbl[i, names(table(as.vector(group.labels[samples])))] <-
                table(as.vector(group.labels[samples]))
            }

            tbl <- tbl / max(tbl)

            stars(tbl, locations=layout, xlim=c(-1.1,1.1), ylim=c(-1.1,1.1),
                  len=0.1, draw.segments=T, col.segments=unique.group.colors,
                  key.loc=c(1.3,-1.1), scale=F)

            plot(stg, layout=layout, vertex.size=0, vertex.label="",
                 xlim=c(-1.1,1.1), ylim=c(-1.1,1.1))

            text(-1.65, 1.1, "Distance Spanning Tree", adj=0, cex=2)
            text(-1.65, 1, "relative to category size", adj=0)

            par(new=T)

            tbl <- matrix(0, ncol(pattern.adj.matrix), length(unique(group.labels)))
            colnames(tbl) <- unique(group.labels)

            for (i in 1:ncol(pattern.adj.matrix))
            {
              samples <- pattern.sample.order[[i]]
              tt <- table(as.vector(group.labels[samples]))

              for (ii in 1:length(tt))
              {
                tt[ii] <- tt[ii] / table(as.vector(group.labels))[names(tt[ii])]
              }

              tbl[i, names(tt)] <- tt
            }

            stars(tbl, locations=layout, xlim=c(-1.1,1.1), ylim=c(-1.1,1.1),
                  len=0.1, draw.segments=T, col.segments=unique.group.colors,
                  key.loc=c(1.3,-1.1), scale=F)
          }
        }
      }
    }
  }

  dirname <- file.path(paste(files.name, "- Results"), "3rd lvl Spot Analysis")

  filename <-
    if (GS.infos.overexpression$filtered)
    {
      file.path(dirname, "Overexpression Networks filtered.pdf")
    } else
    {
      file.path(dirname, "Overexpression Networks.pdf")
    }

  util.info("Writing:", filename)
  pdf(filename, 29.7/2.54, 21/2.54)
  plot.set.list.networks(set.list=GS.infos.overexpression, main="Sample-Overexpression")
  dev.off()

  filename <-
    if (GS.infos.underexpression$filtered)
    {
      file.path(dirname, "Underexpression Networks filtered.pdf")
    } else
    {
      file.path(dirname, "Underexpression Networks.pdf")
    }

  util.info("Writing:", filename)
  pdf(filename, 29.7/2.54, 21/2.54)
  plot.set.list.networks(set.list=GS.infos.underexpression, main="Sample-Underexpression")
  dev.off()

  filename <- file.path(dirname, "K-Means Cluster Networks.pdf")
  util.info("Writing:", filename)
  pdf(filename, 29.7/2.54, 21/2.54)
  plot.set.list.networks(set.list=GS.infos.kmeans, main="K-Means Cluster")
  dev.off()

  filename <- file.path(dirname, "Group Overexpression Networks.pdf")
  util.info("Writing:", filename)
  pdf(filename, 29.7/2.54, 21/2.54)
  plot.set.list.networks(set.list=GS.infos.group.overexpression, main="Group Overexpression")
  dev.off()
}
