pipeline.3rdLvlSummarySheets <- function()
{
  plot.set.list.reports <- function(set.list, main)
  {
    N.spots <- length(set.list$spots)

    ### Binary and Homogeneity maps for groups in the spots
    thresh.global <- sd(as.vector(set.list$spotdata))

    B.spot.group.map.global <- matrix(0,
                                      length(set.list$spots),
                                      length(unique(group.labels)),
                                      dimnames=list(names(set.list$spots),
                                                    unique(group.labels)))

    H.spot.group.map <- matrix(NA,
                               length(set.list$spots),
                               length(unique(group.labels)),
                               dimnames=list(names(set.list$spots),
                                             unique(group.labels)))

    S.spot.group.map <- matrix(NA,
                               length(set.list$spots),
                               length(unique(group.labels)),
                               dimnames=list(names(set.list$spots),
                                             unique(group.labels)))

    for (s in names(set.list$spots))
    {
      low.level <- -thresh.global
      high.level <- thresh.global

      for (g in seq_along(unique(group.labels)))
      {
        gr <- unique(group.labels)[g]
        gr.mem <- which(group.labels == gr)

        p <- hist(set.list$spotdata[s,gr.mem],
                  breaks=c(min(set.list$spotdata), low.level, high.level,
                           max(set.list$spotdata)), plot=FALSE)$counts / length(gr.mem)
        p <- p[p != 0]

        H <- - sum(p * log2(p))
        H <- H / log2(3)

        H.spot.group.map[s, gr] <- H
      }

      B.spot.group.map.global[s,which(tapply(set.list$spotdata[s,], group.labels, mean)[unique(group.labels)] > thresh.global)] <- 1
      B.spot.group.map.global[s, which(tapply(set.list$spotdata[s,], group.labels, mean)[unique(group.labels)] < -thresh.global)] <- -1

      x <- set.list$spotdata[s,]
      x <- x-min(x)
      p <- x / sum(x)
      p <- p[p!=0]

      H <- sum(-p * log2(p))
      H <- H / log2(length(unique(group.labels)))

      S.spot.group.map[s,] <- B.spot.group.map.global[s,] * H
    }

    ### Binary matrix for spot-to-sample assignment
    sample.spot.matrix <- matrix(0, length(set.list$spots), ncol(metadata))
    rownames(sample.spot.matrix) <- names(set.list$spots)

    for (s in names(set.list$spots))
    {
      if (main %in% c("Sample-Underexpression"))
      {
        sample.spot.matrix[s,] <- as.numeric(set.list$spotdata[s,] < -thresh.global)
      } else
      {
        sample.spot.matrix[s,] <- as.numeric(set.list$spotdata[s,] > thresh.global)
      }
    }

    ## plot overview
    par(mfrow=c(1,2))
    par(mar=c(8, 3, 6, 3))

    col <- if(main!="D-Clusters") color.palette.portraits(1000) else colorRampPalette(c("blue2","white","red2"))(1000)
    image(matrix(set.list$overview.map, preferences$dim.1stLvlSom, preferences$dim.1stLvlSom),
          axes=FALSE, col=col, main=paste("Overview map,",main), cex.main=1.5)

    box()

    image(matrix(set.list$overview.mask, preferences$dim.1stLvlSom, preferences$dim.1stLvlSom),
          axes=FALSE, col = "darkgreen", main="Detected spots", cex.main=1)

    box()

    par(new=TRUE)

    plot(0, type="n", axes=TRUE, xlab="", ylab="", xlim=c(0,preferences$dim.1stLvlSom),
         ylim=c(0,preferences$dim.1stLvlSom), xaxs="i", yaxs="i", las=1)

    points(do.call(rbind, lapply(set.list$spots, function(x) x$position)), pch=16, cex=2.6, col="black")
    points(do.call(rbind, lapply(set.list$spots, function(x) x$position)), pch=1, cex=2.6, col="white")
    text(do.call(rbind, lapply(set.list$spots, function(x) x$position)), names(set.list$spots), col="white", cex=0.8)

    ## plot spot profiles
    layout(matrix(1:18, 6, 3, byrow=TRUE), widths=c(1,4,1.5))

    for (i in 1:N.spots)
    {
      samples <- which(sample.spot.matrix[i,] == 1)
      mask <- set.list$spots[[i]]$mask

      par(mar=c(0.5,3,0.5,1))

      image(matrix(mask, preferences$dim.1stLvlSom, preferences$dim.1stLvlSom),
            axes=FALSE, col ="darkgreen")

      axis(2, 0.95, names(set.list$spots[i]), las=2, tick=FALSE, cex.axis=1.6)
      box()

      par(mar=c(0.5,3,0.5,1))

      barplot.x <- barplot(set.list$spotdata[i,], names.arg=NA,
                           col=group.colors, main="", las=2, cex.main=1,
                           cex.lab=1, cex.axis=1, cex.names=0.8,
                           border=if (ncol(indata) < 80) "black" else NA,
                           ylim=range(set.list$spotdata)*c(1,1.4))

      abline(h=thresh.global*c(-1,1), lty=2, col="gray")
      label.string <- as.character(B.spot.group.map.global[i,])
      label.string <- sub("-1", "-", label.string, fixed=TRUE)
      label.string <- sub("1", "+", label.string)
      label.string <- sub("0", ".", label.string)
      label.x <- tapply(barplot.x, group.labels, mean)[unique(group.labels)]
      text(label.x, max(set.list$spotdata)*1.2, label.string, col=groupwise.group.colors,cex=2.5)

      par(mar=c(0.5,0,0.5,0))
      plot(0, type="n", axes=FALSE, xlab="", ylab="", xlim=c(0,1), ylim=c(0,1), xaxs="i", yaxs="i")
      text(0, 0.9, paste("# samples expressing :", length(samples)), adj=0)

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
    ylim <- range( unlist( by( t(set.list$spotdata), group.labels, range )[unique(group.labels)] ) )
   
    for (i in 1:N.spots)
    {
      par(mar=c(2,3,0.5,1))
      
      image(matrix(set.list$spots[[i]]$mask, preferences$dim.1stLvlSom, preferences$dim.1stLvlSom),
            axes=FALSE, col ="darkgreen")
      
      axis(2, 0.95, names(set.list$spots[i]), las=2, tick=FALSE, cex.axis=1.6)
      box()
      
      par(mar=c(2,3,0.5,1))
      
      boxplot(  tapply( set.list$spotdata[i,], group.labels, c )[unique(group.labels)],
                col=groupwise.group.colors, las=1, ylim=ylim )
      
      abline(h=thresh.global*c(-1,0,1), lty=2, col="gray")

      label.string <- as.character(B.spot.group.map.global[i,])
      label.string <- sub("-1", "-", label.string, fixed=TRUE)
      label.string <- sub("1", "+", label.string)
      label.string <- sub("0", ".", label.string)

      text( 1:length(unique(group.labels)), par("usr")[4]-(par("usr")[4]-par("usr")[3])*0.1, label.string, col=groupwise.group.colors,cex=3)
    }
    
    ### Spot homogeneity map
    map <- B.spot.group.map.global*(1-H.spot.group.map)
    map <- map[nrow(map):1,,drop=FALSE]

    layout(matrix(c(1, 2), 1, 2), c(3, 1), 1)
    par(mar=c(5,4,4,0))

    image(t(map), col=colorRampPalette(c("darkblue","white","darkred"))(1000),
          main="Spot homogeneity map", zlim=c(-1,1), axes=FALSE, cex.main=1.5)

    box()
    axis(1, seq(0,1,length.out=ncol(map)), colnames(map), las=2, cex.axis=1.1)
    axis(2, seq(0,1,length.out=nrow(map)), rownames(map), las=2, cex.axis=1.1)

    par(mar=c(0,0,4,0))
    plot(0,type="n",axes=FALSE,xlab="",ylab="")

    legend("top", c("Overexpressed, homogenous","Overexpressed, heterogenous",
                    "Not differentially expressed","Underexpressed, heterogenous",
                    "Underexpressed, homogenous"),
           pch=c(15,15,0,15,15), col=c("darkred","pink","black","lightblue","darkblue"))

    map = 1-H.spot.group.map
    map[which(B.spot.group.map.global<=0)] = 0

    layout(matrix(c(1, 2), 1, 2), c(2, 1), 1)
    par(mar=c(5, 4, 4, 2))

    image(matrix(set.list$overview.mask, preferences$dim.1stLvlSom), col="mistyrose", axes=FALSE, main="Spot homogeneity map - overexpressed spots", cex.main=1.5)
    par(new=TRUE)

    plot(0, type="n", xlab="", ylab="", axes=FALSE, xlim=c(0,preferences$dim.1stLvlSom),
         ylim=c(0,preferences$dim.1stLvlSom), xaxs="i", yaxs="i")

    box()

    stars(map, locations=t(sapply(set.list$spots, function(x) x$position)), len=2,
          draw.segments=TRUE,  scale=FALSE, add=TRUE, col.segments=groupwise.group.colors, labels=NULL)

    text(sapply(set.list$spots, function(x) x$position)[1,],
         sapply(set.list$spots, function(x) x$position)[2,]+preferences$dim.1stLvlSom*0.05,
         names(set.list$spots), col="gray50")

    par(mar=c(0,0,4,0))
    plot(0,type="n",axes=FALSE,xlab="",ylab="")

    legend("topleft", legend=unique(group.labels), cex=1, col=groupwise.group.colors,
           pch=15, pt.cex=1.4)   
    
    map = 1-H.spot.group.map
    map[which(B.spot.group.map.global>=0)] <- 0

    layout(matrix(c(1, 2), 1, 2), c(2, 1), 1)
    par(mar=c(5, 4, 4, 2))

    image(matrix(set.list$overview.mask, preferences$dim.1stLvlSom), col="lightsteelblue1", main="Spot homogeneity map - underexpressed spots", cex.main=1.5)
    par(new=TRUE)

    plot(0, type="n", xlab="", ylab="", axes=FALSE, xlim=c(0,preferences$dim.1stLvlSom),
         ylim=c(0,preferences$dim.1stLvlSom), xaxs="i", yaxs="i") #, main="Group Association", cex.main=2.5

    box()

    stars(map, locations=t(sapply(set.list$spots, function(x) x$position)),
          len=2, draw.segments=TRUE,  scale=FALSE, add=TRUE, col.segments=groupwise.group.colors,
          labels=NULL)

    text(sapply(set.list$spots, function(x) x$position)[1,],
         sapply(set.list$spots, function(x) x$position)[2,]+preferences$dim.1stLvlSom*0.05,
         names(set.list$spots), col="gray50")

    par(mar=c(0,0,4,0))
    plot(0,type="n",axes=FALSE,xlab="",ylab="")

    legend("topleft", legend=unique(group.labels), cex=1, col=groupwise.group.colors,
           pch=15, pt.cex=1.4)   
  
    ### Group association barplots
    if (length(unique(group.labels)) > 1)
    {
      spot.group.assoc <- apply(sample.spot.matrix, 1, function(x)
      {
        tapply(x, group.labels, sum) / table(group.labels)
      })

      spot.group.assoc <- spot.group.assoc[unique(group.labels) , ,drop=FALSE]
      colnames(spot.group.assoc) <- names(set.list$spots)

      
      spot.goup.mean.number <- tapply(colSums(sample.spot.matrix), group.labels, mean)[unique(group.labels)]

      layout(matrix(c(1, 2), 1, 2), c(2, 1), 1)
      par(mar=c(5, 4, 4, 2))

      barplot(spot.group.assoc[,ncol(spot.group.assoc):1],
              main="Association of the groups to the spots", cex.main=1.5, cex.axis=2, cex.names=1,
              col=groupwise.group.colors, horiz=TRUE, las=1)

      par(mar=c(5, 1, 4, 2))

      plot(0, main="", cex.main=1, type="n", axes=FALSE, xlab="",
           ylab="", xlim=c(0,1), ylim=c(0,1), xaxs="i", yaxs="i")

      legend("topleft",  legend=paste(unique(group.labels), "(", round(spot.goup.mean.number,1), ")"),
             col=groupwise.group.colors, pch=15, pt.cex=1.4, title="< #spots >")

      layout(matrix(c(1, 2), 1, 2), c(2, 1), 1)
      par(mar=c(5, 4, 4, 2))

      image(matrix(set.list$overview.mask, preferences$dim.1stLvlSom), col="gray90", axes=FALSE,
            main="Association of the groups to the spots", cex.main=1.5)
      par(new=TRUE)
      plot(0, type="n", xlab="", ylab="", axes=FALSE, xlim=c(0,preferences$dim.1stLvlSom),
           ylim=c(0,preferences$dim.1stLvlSom), xaxs="i", yaxs="i")

      box()

      for (i in seq_along(set.list$spots))
      {
        stars(t(spot.group.assoc[,i]), locations=set.list$spots[[i]]$position,
              len=2, draw.segments=TRUE,  scale=FALSE, add=TRUE, col.segments=groupwise.group.colors)

        par(fg="gray", lty=2)

        stars(t(c(0.5)), locations=set.list$spots[[i]]$position, len=2,
              draw.segments=TRUE,  scale=FALSE, add=TRUE, col.segments=NA)

        par(fg="black", lty="solid")

        text(set.list$spots[[i]]$position[1], set.list$spots[[i]]$position[2]+preferences$dim.1stLvlSom*0.05,
             names(set.list$spots)[i], col="gray50")
      }

      par(mar=c(5, 1, 4, 2))

      plot(0, type="n", axes=FALSE, xlab="", ylab="", xlim=c(0,1), ylim=c(0,1),
           xaxs="i", yaxs="i")

      legend("topleft", legend=unique(group.labels), col=groupwise.group.colors, pch=15, pt.cex=1.4)
    }
  }

  dirname <- file.path(paste(files.name, "- Results"), "3rd lvl Spot Analysis")

  filename <-
    if (spot.list.overexpression$filtered)
    {
      file.path(dirname, "Spot Report - Overexpression Spots filtered.pdf")
    } else
    {
      file.path(dirname, "Spot Report - Overexpression Spots.pdf")
    }

  util.info("Writing:", filename)
  pdf(filename, 29.7/2.54, 21/2.54)
  plot.set.list.reports(set.list=spot.list.overexpression, main="Overexpression Spots")
  dev.off()

  filename <-
    if (spot.list.underexpression$filtered)
    {
      file.path(dirname, "Spot Report - Underexpression Spots filtered.pdf")
    } else
    {
      file.path(dirname, "Spot Report - Underexpression Spots.pdf")
    }

  util.info("Writing:", filename)
  pdf(filename, 29.7/2.54, 21/2.54)
  plot.set.list.reports(set.list=spot.list.underexpression, main="Underexpression Spots")
  dev.off()

  filename <- file.path(dirname, "Spot Report - K-Means Clusters.pdf")
  util.info("Writing:", filename)
  pdf(filename, 29.7/2.54, 21/2.54)
  plot.set.list.reports(set.list=spot.list.kmeans, main="K-Means Clusters")
  dev.off()

  if (length(unique(group.labels)) > 1)
  {
    filename <- file.path(dirname, "Spot Report - Group Overexpression Spots.pdf")
    util.info("Writing:", filename)
    pdf(filename, 29.7/2.54, 21/2.54)
    plot.set.list.reports(set.list=spot.list.group.overexpression, main="Group Overexpression Spots")
    dev.off()
  }
  
  filename <- file.path(dirname, "Spot Report - D-Clusters.pdf")
  util.info("Writing:", filename)
  pdf(filename, 29.7/2.54, 21/2.54)
  plot.set.list.reports(set.list=spot.list.dmap, main="D-Clusters")
  dev.off()
}
