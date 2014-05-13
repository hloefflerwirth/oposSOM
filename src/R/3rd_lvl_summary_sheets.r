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

      for (g in 1:length(unique(group.labels)))
      {
        gr <- unique(group.labels)[g]
        gr.mem <- which(group.labels == gr)

        p <- hist(set.list$spotdata[s,gr.mem],
                  breaks=c(min(set.list$spotdata), low.level, high.level,
                           max(set.list$spotdata)), plot=F)$counts / length(gr.mem)
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

    image(matrix(set.list$overview.map, preferences$dim.1stLvlSom, preferences$dim.1stLvlSom),
          axes=F, col = colramp(1000), main=main, cex.main=1.5)

    box()

    image(matrix(set.list$overview.mask, preferences$dim.1stLvlSom, preferences$dim.1stLvlSom),
          axes=F, col = "darkgreen", main="Spots", cex.main=1.5)

    box()

    par(new=T)

    plot(0, type="n", axes=T, xlab="", ylab="", xlim=c(0,preferences$dim.1stLvlSom),
         ylim=c(0,preferences$dim.1stLvlSom), xaxs="i", yaxs="i", las=1)

    points(do.call(rbind, lapply(set.list$spots, function(x) x$position)), pch=16, cex=3, col="black")
    points(do.call(rbind, lapply(set.list$spots, function(x) x$position)), pch=1, cex=3, col="white")
    text(do.call(rbind, lapply(set.list$spots, function(x) x$position)), names(set.list$spots), col="white")

    ## plot single spot overview
    layout(matrix(1:18, 6, 3, byrow=T), widths=c(1,4,1.5))

    for (i in 1:N.spots)
    {
      samples <- which(sample.spot.matrix[i,] == 1)
      mask <- set.list$spots[[i]]$mask

      par(mar=c(0.5,3,0.5,1))

      image(matrix(mask, preferences$dim.1stLvlSom, preferences$dim.1stLvlSom),
            axes=F, col ="darkgreen")

      axis(2, 0.95, names(set.list$spots[i]), las=2, tick=F, cex.axis=1.6)
      box()

      par(mar=c(0.5,3,0.5,1))

      barplot.x <- barplot(set.list$spotdata[i,], names.arg=NA,
                           col=group.colors, main="", las=2, cex.main=1,
                           cex.lab=1, cex.axis=1, cex.names=0.8,
                           border=if (ncol(indata) < 80) "black" else NA,
                           ylim=range(set.list$spotdata)*c(1,1.4))

      abline(h=thresh.global, lty=2, col="gray")
      abline(h=-thresh.global, lty=2, col="gray")
      label.string = as.character(B.spot.group.map.global[i,])
      label.string = sub("-1", "-", label.string, fixed=T)
      label.string = sub("1", "+", label.string)
      label.string = sub("0", ".", label.string)
      label.x = tapply(barplot.x, group.labels, mean)[unique(group.labels)]
      text(label.x, max(set.list$spotdata)*1.2, label.string, col=groupwise.group.colors,cex=2.5)

      par(mar=c(0.5,0,0.5,0))
      plot(0, type="n", axes=F, xlab="", ylab="", xlim=c(0,1), ylim=c(0,1), xaxs="i", yaxs="i")
      text(0, 0.9, paste("# samples expressing :", length(samples)), adj=0)

      if (length(samples) > 0)
      {
        t <- table(as.vector(group.labels[samples]))[unique(group.labels)]
        t <- t[!is.na(t)]

        for (ii in 1:length(t))
        {
          t[ii] <- t[ii] / table(as.vector(group.labels))[names(t[ii])]
        }

        t <- round(100 * t)

        t2 <- table(as.vector(group.labels[samples]))[unique(group.labels)]
        t2 <- t2[!is.na(t2)]

        for (ii in 1:length(t))
        {
          text(0, 0.8  - ii*0.1,
               paste(names(t)[ii], ": # =", t2[ii], " -> ", t[ii], "%"),
               adj=0, cex=1, col=group.colors[which(group.labels == names(t)[ii])[1]])
        }
      }
    }

    ### Binary spot-group-heatmaps
    map <- B.spot.group.map.global*(1-H.spot.group.map)
    map <- map[nrow(map):1,,drop=F]

    layout(matrix(c(1, 2), 1, 2), c(3, 1), 1)
    par(mar=c(5,4,4,0))

    image(t(map), col=colorRampPalette(c("darkblue","white","darkred"))(1000),
          main="Spot homogeneity map", zlim=c(-1,1), axes=F, cex.main=2.5)

    box()
    axis(1, seq(0,1,length.out=ncol(map)), colnames(map), las=2, cex.axis=1.1)
    axis(2, seq(0,1,length.out=nrow(map)), rownames(map), las=2, cex.axis=1.1)

    par(mar=c(0,0,4,0))
    plot(0,type="n",axes=F,xlab="",ylab="")

    legend("top", c("Overexpressed, homogenous","Overexpressed, heterogenous",
                    "Not differentially expressed","Underexpressed, heterogenous",
                    "Underexpressed, homogenous"),
           pch=c(15,15,0,15,15), col=c("darkred","pink","black","lightblue","darkblue"))

    map = 1-H.spot.group.map
    map[which(B.spot.group.map.global<=0)] = 0

    layout(matrix(c(1, 2), 1, 2), c(2, 1), 1)
    par(mar=c(5, 4, 4, 2))

    image(matrix(set.list$overview.mask, preferences$dim.1stLvlSom), col="mistyrose", axes=F)
    par(new=T)

    plot(0, type="n", xlab="", ylab="", axes=F, xlim=c(0,preferences$dim.1stLvlSom),
         ylim=c(0,preferences$dim.1stLvlSom), xaxs="i", yaxs="i") #, main="Group Association", cex.main=2.5

    box()

    stars(map, locations=t(sapply(set.list$spots, function(x) x$position)), len=2,
          draw.segments=T,  scale=F, add=T, col.segments=groupwise.group.colors, labels=NULL)

    text(sapply(set.list$spots, function(x) x$position)[1,],
         sapply(set.list$spots, function(x) x$position)[2,]+preferences$dim.1stLvlSom*0.05,
         names(set.list$spots), col="gray50")

    par(mar=c(0,0,4,0))
    plot(0,type="n",axes=F,xlab="",ylab="")

    legend("topleft", legend=unique(group.labels), cex=1.3, col=groupwise.group.colors,
           pch=15, pt.cex=2, bty="n")

    map = 1-H.spot.group.map
    map[which(B.spot.group.map.global>=0)] <- 0

    layout(matrix(c(1, 2), 1, 2), c(2, 1), 1)
    par(mar=c(5, 4, 4, 2))

    image(matrix(set.list$overview.mask, preferences$dim.1stLvlSom), col="lightsteelblue1", axes=F)
    par(new=T)

    plot(0, type="n", xlab="", ylab="", axes=F, xlim=c(0,preferences$dim.1stLvlSom),
         ylim=c(0,preferences$dim.1stLvlSom), xaxs="i", yaxs="i") #, main="Group Association", cex.main=2.5

    box()

    stars(map, locations=t(sapply(set.list$spots, function(x) x$position)),
          len=2, draw.segments=T,  scale=F, add=T, col.segments=groupwise.group.colors,
          labels=NULL)

    text(sapply(set.list$spots, function(x) x$position)[1,],
         sapply(set.list$spots, function(x) x$position)[2,]+preferences$dim.1stLvlSom*0.05,
         names(set.list$spots), col="gray50")

    par(mar=c(0,0,4,0))
    plot(0,type="n",axes=F,xlab="",ylab="")

    legend("topleft", legend=unique(group.labels), cex=1.3, col=groupwise.group.colors,
           pch=15, pt.cex=2, bty="n")

    ### Group association barplots
    if (length(unique(group.labels)) > 1)
    {
      spot.group.assoc <- apply(sample.spot.matrix, 1, function(x)
      {
        tapply(x, group.labels, sum) / table(group.labels)
      })

      spot.group.assoc <- spot.group.assoc[unique(group.labels) , ,drop=F]
      colnames(spot.group.assoc) <- names(set.list$spots)

      spot.goup.mean.number <- sapply(spot.list.samples, function(x)
      {
        sum (sapply(x$spots, function(y)
        {
          if (main %in% c("Sample-Underexpression"))
          {
            if (y$type == "underexpressed") 1 else 0
          }  else
          {
            if (y$type == "overexpressed") 1 else 0
          }
        }))
      })

      spot.goup.mean.number <- tapply(spot.goup.mean.number, group.labels, mean)[unique(group.labels)]

      layout(matrix(c(1, 2), 1, 2), c(2, 1), 1)
      par(mar=c(5, 4, 4, 2))

      barplot(spot.group.assoc[,ncol(spot.group.assoc):1],
              main="Group Association", cex.main=2.5, cex.axis=2, cex.names=1,
              col=groupwise.group.colors, horiz=T, las=1)

      par(mar=c(5, 1, 4, 2))

      plot(0, main="< #spots >", cex.main=2.5, type="n", axes=F, xlab="",
           ylab="", xlim=c(0,1), ylim=c(0,1), xaxs="i", yaxs="i")

      legend("topleft",
             legend=paste(unique(group.labels), "(", round(spot.goup.mean.number,1), ")"),
             cex=1.3, col=groupwise.group.colors, pch=15, pt.cex=2, bty="n")

      layout(matrix(c(1, 2), 1, 2), c(2, 1), 1)
      par(mar=c(5, 4, 4, 2))

      image(matrix(set.list$overview.mask, preferences$dim.1stLvlSom), col="gray90", axes=F)
      par(new=T)
      plot(0, type="n", xlab="", ylab="", axes=F, xlim=c(0,preferences$dim.1stLvlSom),
           ylim=c(0,preferences$dim.1stLvlSom), xaxs="i", yaxs="i") #, main="Group Association", cex.main=2.5

      box()

      for (i in 1:length(set.list$spots))
      {
        stars(t(spot.group.assoc[,i]), locations=set.list$spots[[i]]$position,
              len=2, draw.segments=T,  scale=F, add=T, col.segments=groupwise.group.colors)

        par(fg="gray", lty=2)

        stars(t(c(0.5)), locations=set.list$spots[[i]]$position, len=2,
              draw.segments=T,  scale=F, add=T, col.segments=NA)

        par(fg="black", lty="solid")

        text(set.list$spots[[i]]$position[1], set.list$spots[[i]]$position[2]+preferences$dim.1stLvlSom*0.05,
             names(set.list$spots)[i], col="gray50")
      }

      par(mar=c(5, 1, 4, 2))

      plot(0, type="n", axes=F, xlab="", ylab="", xlim=c(0,1), ylim=c(0,1),
           xaxs="i", yaxs="i")

      legend("topleft", legend=unique(group.labels), cex=1.3, col=groupwise.group.colors,
             pch=15, pt.cex=2, bty="n")

      group.spot.mean.number <- apply(ceiling(spot.group.assoc), 2, sum)

      layout(matrix(c(1, 2), 1, 2), c(2, 1), 1)
      par(mar=c(5, 8, 4, 2))

      barplot(t(spot.group.assoc)[,nrow(spot.group.assoc):1],
              main="Spot Association", cex.main=2.5, cex.axis=2, cex.names=0.8,
              col=colramp(N.spots), horiz=T, las=1)

      par(mar=c(5, 1, 4, 2))

      plot(0, main="#groups", cex.main=2.5, type="n", axes=F, xlab="", ylab="",
           xlim=c(0,1), ylim=c(0,1), xaxs="i", yaxs="i")

      legend("topleft",
             legend=paste(names(set.list$spots), "(", group.spot.mean.number, ")"),
             cex=1.3, col=colramp(N.spots), pch=15, pt.cex=2, bty="n")
    }

    if (preferences$geneset.analysis)
    {
      sig.gs.types <- matrix(0,
                             length(set.list$spots),
                             length(unique(gs.def.list.categories)),
                             dimnames=list(names(set.list$spots),
                                           unique(gs.def.list.categories)))

      for (i in 1:length(set.list$spots))
      {
        x <- set.list$spots[[i]]

        sig.gs <- names(which(x$Fisher.p < 1e-4))
        sig.gs.types.table <- table(gs.def.list.categories[sig.gs])

        sig.gs.types[i, names(sig.gs.types.table)] <- sig.gs.types.table
      }

      for (i in 1:ncol(sig.gs.types))
      {
        sig.gs.types[,i] <-
          sig.gs.types[,i] / sum(gs.def.list.categories == colnames(sig.gs.types)[i]) * 100
      }

      layout(matrix(c(1, 2), 1, 2), c(2, 1), 1)
      par(mar=c(5, 4, 4, 2))

      barplot(t(sig.gs.types)[,nrow(sig.gs.types):1], main="%Significant Genesets",
              cex.main=2.5, cex.axis=2, cex.names=1,
              col=colramp(length(unique(gs.def.list.categories))),
              horiz=T, las=1)

      par(mar=c(5, 1, 4, 2))

      plot(0, main="", cex.main=2.5, type="n", axes=F, xlab="", ylab="",
           xlim=c(0,1), ylim=c(0,1), xaxs="i", yaxs="i")

      legend("topleft", legend=unique(gs.def.list.categories), cex=1.3,
             col=colramp(length(unique(gs.def.list.categories))), pch=15,
             pt.cex=2, bty="n")

      sig.gs.types = sig.gs.types / max(sig.gs.types)

      layout(matrix(c(1, 2), 1, 2), c(2, 1), 1)
      par(mar=c(5, 4, 4, 2))

      image(matrix(set.list$overview.mask, preferences$dim.1stLvlSom), col="gray90", axes=F)
      par(new=T)

      plot(0, type="n", xlab="", ylab="", axes=F, xlim=c(0,preferences$dim.1stLvlSom),
           ylim=c(0,preferences$dim.1stLvlSom), xaxs="i", yaxs="i") #, main="Group Association", cex.main=2.5

      box()

      for (i in 1:length(set.list$spots))
      {
        stars(t(sig.gs.types[i,]), locations=set.list$spots[[i]]$position, len=2,
              draw.segments=T,  scale=F, add=T,
              col.segments=colramp(length(unique(gs.def.list.categories))))

        par(fg="gray", lty=2)

        stars(t(c(0.5)), locations=set.list$spots[[i]]$position, len=2,
              draw.segments=T,  scale=F, add=T, col.segments=NA)

        par(fg="black", lty="solid")

        text(set.list$spots[[i]]$position[1],
             set.list$spots[[i]]$position[2]+preferences$dim.1stLvlSom*0.05,
             names(set.list$spots)[i], col="gray50")
      }

      par(mar=c(5, 1, 4, 2))

      plot(0, type="n", axes=F, xlab="", ylab="", xlim=c(0,1), ylim=c(0,1),
           xaxs="i", yaxs="i")

      legend("topleft",
             legend=unique(gs.def.list.categories),
             cex=1.3, col=colramp(length(unique(gs.def.list.categories))),
             pch=15, pt.cex=2, bty="n")

    }
  }

  dirname <- file.path(paste(files.name, "- Results"), "3rd lvl Spot Analysis")

  filename <-
    if (spot.list.overexpression$filtered)
    {
      file.path(dirname, "Overexpression Spot Report filtered.pdf")
    } else
    {
      file.path(dirname, "Overexpression Spot Report.pdf")
    }

  util.info("Writing:", filename)
  pdf(filename, 29.7/2.54, 21/2.54)
  plot.set.list.reports(set.list=spot.list.overexpression, main="Sample-Overexpression")
  dev.off()

  filename <-
    if (spot.list.underexpression$filtered)
    {
      file.path(dirname, "Underexpression Spot Report filtered.pdf")
    } else
    {
      file.path(dirname, "Underexpression Spot Report.pdf")
    }

  util.info("Writing:", filename)
  pdf(filename, 29.7/2.54, 21/2.54)
  plot.set.list.reports(set.list=spot.list.underexpression, main="Sample-Underexpression")
  dev.off()

  filename <- file.path(dirname, "K-Means Cluster Report.pdf")
  util.info("Writing:", filename)
  pdf(filename, 29.7/2.54, 21/2.54)
  plot.set.list.reports(set.list=spot.list.kmeans, main="K-Means Cluster")
  dev.off()

  if (length(unique(group.labels)) > 1)
  {
    filename <- file.path(dirname, "Group Overexpression Report.pdf")
    util.info("Writing:", filename)
    pdf(filename, 29.7/2.54, 21/2.54)
    plot.set.list.reports(set.list=spot.list.group.overexpression, main="Group Overexpression")
    dev.off()
  }
}
