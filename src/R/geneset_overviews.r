pipeline.genesetOverviews <- function()
{
  filename <- file.path(paste(files.name, "- Results"), "Geneset Analysis", "0verview Heatmaps.pdf")
  util.info("Writing:", filename)
  pdf(filename, 21/2.54, 21/2.54)

  ### GSZ heatmaps
  spots.GSZ.scores <- do.call(cbind, lapply(spot.list.samples, function(x)
  {
    apply(do.call(cbind,  lapply(x$spots, function(y) { y$GSZ.score[names(gs.def.list)] })),
          1, function(z) { if (max(z) > -min(z)) { max(z) } else { min(z) } })
  }))

  summary.spots.fisher.p <-
    sapply(spot.list.overexpression$spots, function(x) { x$Fisher.p[names(gs.def.list)] })


  for (i in names(table(gs.def.list.categories)))
  {
    category.GST.scores <- spots.GSZ.scores[which(gs.def.list.categories == i) , ,drop=F]

    if (nrow(category.GST.scores) > 60)
    {
      category.GST.order <- apply(category.GST.scores, 2, order, decreasing=TRUE)
      category.GST.order <- as.vector(t(category.GST.order))

      j <- 1

      while (length(unique(category.GST.order[1:(ncol(indata)*j)])) < 60)
      {
        j <- j + 1
      }

      top.scores <- category.GST.scores[unique(category.GST.order[1:(ncol(indata)*j)]),]
    } else
    {
      top.scores <- category.GST.scores
    }

    top.scores[which(top.scores > quantile(top.scores, 0.99))] <- quantile(top.scores, 0.99)

    if (nrow(top.scores) == 1)
    {
      # cracky workaround: heatmap demands 2 rows
      top.scores <- rbind(top.scores, top.scores)
    }

    rownames(top.scores) <-
      paste(LETTERS[apply(summary.spots.fisher.p[rownames(top.scores), ,drop=F], 1, which.min)],
            rownames(top.scores), sep="  ")

    colnames(top.scores) <- colnames(indata)

    heatmap.wrap(x=top.scores, cex.main=2,
                 col=colorRampPalette(c("blue4","blue","gray90","orange","red4"))(1000),
                 mar=c(10,20), scale="n", zlim=max(max(top.scores),-min(top.scores))*c(-1,1),
                 ColSideColors=group.colors, cexDend=0.6)

    par(new=TRUE, mar=c(3.5,29,35.5,2))

    image(matrix(c(1:1000), 1000, 1), axes=FALSE,
          col=colorRampPalette(c("blue4","blue","gray90","orange","red4"))(1000))

    box()
    axis(1, c(0,0.5,1), round(max(max(top.scores),-min(top.scores))*c(-1,0,1)), cex.axis=1.4)
    mtext("GSZ score", cex=2, line=32)
    mtext(i, cex=1.2, line=30)

    heatmap.wrap(x=top.scores, cex.main=2,
                 col=colorRampPalette(c("blue4","blue","gray90","orange","red4"))(1000),
                 mar=c(10,20), scale="n",  zlim=max(max(top.scores),-min(top.scores))*c(-1,1),
                 ColSideColors=group.colors, cexDend=0.6, Colv=NA)

    par(new=TRUE, mar=c(3.5,29,35.5,2))

    image(matrix(c(1:1000), 1000, 1), axes=FALSE,
          col=colorRampPalette(c("blue4","blue","gray90","orange","red4"))(1000))

    box()
    axis(1, c(0,0.5,1), round(max(max(top.scores),-min(top.scores))*c(-1,0,1)), cex.axis=1.4)
    mtext("GSZ score", cex=2, line=33.5)
    mtext(i, cex=1.2, line=32)

    o <- unlist(sapply(unique(group.labels), function(gr)
    {
      idx <- names(group.labels)[which(group.labels == gr)]

      if (length(idx) > 1)
      {
        hc <- hclust(dist(t(top.scores[,idx])))
        return(hc$labels[hc$order])
      }
      return(idx)
    }))

    heatmap.wrap(x=top.scores[,o], cex.main=2,
                 col=colorRampPalette(c("blue4","blue","gray90","orange","red4"))(1000),
                 mar=c(10,20), scale="n",
                 zlim=max(max(top.scores),-min(top.scores))*c(-1,1),
                 ColSideColors=group.colors[o], cexDend=0.6, Colv=NA)

    par(new=TRUE, mar=c(3.5,29,35.5,2) )

    image(matrix(c(1:1000), 1000, 1), axes=FALSE,
          col=colorRampPalette(c("blue4","blue","gray90","orange","red4"))(1000))

    box()
    axis( 1, c(0,0.5,1), round(max(max(top.scores),-min(top.scores))*c(-1,0,1)), cex.axis=1.4 )
    mtext( "GSZ score", cex=2, line=33.5 )
    mtext( i, cex=1.2, line=32 )
  }

  ### Overrepresentation heatmaps
  spots.fisher.p <- do.call(cbind, lapply(spot.list.samples, function(x)
  {
    apply(do.call(cbind, lapply(x$spots, function(y)
    {
      if (y$type == "overexpressed")
      {
        y$Fisher.p[names(gs.def.list)]
      }
    })), 1, min)
  }))

  for (i in names(table(gs.def.list.categories)))
  {
    category.fisher.p <- spots.fisher.p[which(gs.def.list.categories == i),,drop=F]

    if (nrow(category.fisher.p) > 60)
    {
      category.fisher.p.order <- apply(category.fisher.p, 2, order, decreasing=FALSE)
      category.fisher.p.order <- as.vector(t(category.fisher.p.order))
      j <- 1

      while (length(unique(category.fisher.p.order[1:(ncol(indata)*j)])) < 60)
      {
        j <- j + 1
      }

      top.scores <- category.fisher.p[unique(category.fisher.p.order[1:(ncol(indata)*j)]),]
    } else
    {
      top.scores <- category.fisher.p
    }

    top.scores <- -log10(top.scores)

    if (sum(is.infinite(top.scores)))
    {
      top.scores[which(is.infinite(top.scores))] <-
        max(top.scores[-which(is.infinite(top.scores))])
    }

    if (nrow(top.scores) == 1)
    {
      # cracky workaround: heatmap demands 2 rows
      top.scores <- rbind(top.scores, top.scores)
    }

    rownames(top.scores) <-
      paste(LETTERS[apply(summary.spots.fisher.p[rownames(top.scores), ,drop=F], 1, which.min)],
            rownames(top.scores), sep="  ")

    # workaround to heatmap failing when too many equal samples are to be sorted
    if (sum(colSums(top.scores) == 0) < ncol(indata) * 0.1)
    {
      heatmap.wrap(x=top.scores, cex.main=2,
                   col=colorRampPalette(c("gray90","orange","red4"))(1000),
                   mar=c(10,20), scale="n", ColSideColors=group.colors, cexDend=0.6)

      par(new=TRUE, mar=c(3.5,29,35.5,2))

      image(matrix(c(1000:1), 1000, 1), axes=FALSE,
            col=colorRampPalette(c("gray90","orange","red4"))(1000))

      box()
      axis(1, c(0,1), c(-round(max(top.scores)), 0), cex.axis=1.4)
      mtext("log(p.value)", cex=2, line=32)
      mtext(i, cex=1.2, line=30)
    }

    heatmap.wrap(x=top.scores, cex.main=2,
                 col=colorRampPalette(c("gray90","orange","red4"))(1000),
                 mar=c(10,20), scale="n", ColSideColors=group.colors, cexDend=0.6, Colv=NA)

    par(new=TRUE, mar=c(3.5,29,35.5,2))
    image(matrix(c(1000:1), 1000, 1), axes=FALSE,
          col=colorRampPalette(c("gray90","orange","red4"))(1000))

    box();
    axis(1, c(0,1), c(-round(max(top.scores)), 0), cex.axis=1.4)
    mtext("log(p.value)", cex=2, line=33.5)
    mtext(i, cex=1.2, line=32)
  }

  ### p Histogram + FDR
  if (preferences$geneset.analysis.exact)
  {
    all.sets.GSZ.p <- unlist(sapply(spot.list.samples, function(x)
    {
      sapply(x$spots, function(x) { x$GSZ.p.value })
    }))

    all.samples.GSZ.p <- unlist(sapply(spot.list.samples, function(x)
    {
      x$GSZ.p.value
    }))

    p <- c(all.sets.GSZ.p, all.samples.GSZ.p)

    suppressWarnings({ fdrtool.result <- fdrtool(p, statistic="pvalue", verbose=FALSE, plot=FALSE) })
    fdr.spot.list.samples <- fdrtool.result$lfdr
    Fdr.spot.list.samples <- fdrtool.result$qval
    n.0.spot.list.samples <- fdrtool.result$param[1,"eta0"]
    perc.DE.spot.list.samples <- 1 - n.0.spot.list.samples

    par(mar=c(5,6,4,5))

    hist(p, bre=20, freq=FALSE, main="p-values (GSZ)", ylab="", xlab="", las=1,
         cex.main=2.5, cex.lab=2, cex.axis=2)

    box()
    mtext("Density", side=2, line=4, cex=2)
    mtext("FDR", side=4, line=4, cex=2)
    mtext(paste("%DE =", round(perc.DE.spot.list.samples, 2)), line=-1.3, cex=1.6)

    abline(h = n.0.spot.list.samples, col="gray", lwd=2)

    par(new=TRUE)
    plot(0, type="n", xlim=c(0,1), ylim=c(0,1), xlab="", ylab="", axes=FALSE)
    axis(4, seq(0, 1, 0.2), seq(0, 1, 0.2), las=1, cex.axis=2)
    o <- order(p)
    o <- o[round(seq(1, length(o), length.out=1000))]
    lines(p[o], Fdr.spot.list.samples[o], lty=2, lwd=2)
    lines(p[o], fdr.spot.list.samples[o], lty=3, lwd=3)

    legend("topright", c("p", expression(eta[0]), "Fdr", "fdr"),
           col=c("black","gray","black","black"), lty=c(1,1,2,3), lwd=c(2,2,2,3))
  }

  all.sets.fisher.p <- unlist(sapply(spot.list.samples, function(x)
  {
    sapply(x$spots, function(x) { x$Fisher.p })
  }))

  p <- all.sets.fisher.p
  length.out <- length(which(p >= 0.99))
  p <- p[p < 0.99]
  p <- c(p, seq(0,1, length.out=length.out))

  suppressWarnings({ fdrtool.result <- fdrtool(p, statistic="pvalue", verbose=FALSE, plot=FALSE) })
  fdr.spot.list.samples <- fdrtool.result$lfdr
  Fdr.spot.list.samples <- fdrtool.result$qval
  n.0.spot.list.samples <- fdrtool.result$param[1,"eta0"]
  perc.DE.spot.list.samples <- 1 - n.0.spot.list.samples

  par(mar=c(5,6,4,5))

  hist(p, bre=20, freq=FALSE, main="p-value (Fisher)", ylab="", xlab="", las=1,
       cex.main=2.5, cex.lab=2, cex.axis=2)

  box()
  mtext("Density", side=2, line=4, cex=2)
  mtext("FDR", side=4, line=4, cex=2)
  mtext(paste("%DE =", round(perc.DE.spot.list.samples ,2)), line=-1.3, cex=1.6)

  abline(h=n.0.spot.list.samples, col="gray", lwd=2)

  par(new=TRUE)
  plot(0, type="n", xlim=c(0,1), ylim=c(0,1), xlab="", ylab="", axes=FALSE)
  axis(4, seq(0, 1, 0.2), seq(0, 1, 0.2), las=1, cex.axis=2)
  o <- order(p)
  lines(p[o], Fdr.spot.list.samples[o], lty=2, lwd=2)
  lines(p[o], fdr.spot.list.samples[o], lty=3, lwd=3)

  legend("topright", c("p", expression(eta[0]), "Fdr", "fdr"),
         col=c("black","gray","black","black"), lty=c(1,1,2,3), lwd=c(2,2,2,3))

  dev.off()
}
