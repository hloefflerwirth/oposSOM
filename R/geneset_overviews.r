pipeline.genesetOverviews <- function(env)
{
  filename <- file.path("Geneset Analysis", "0verview Heatmaps.pdf")
  util.info("Writing:", filename)
  pdf(filename, 21/2.54, 21/2.54, useDingbats=FALSE)

  ### GSZ heatmaps
  summary.spots.fisher.p <-
    sapply(env[[paste("spot.list.",env$preferences$standard.spot.modules,sep="")]]$spots, function(x) { x$Fisher.p[names(env$gs.def.list)] })


  for (i in names(table(sapply(env$gs.def.list, function(x) { x$Type }))))
  {
    category.GST.scores <- env$samples.GSZ.scores[which(sapply(env$gs.def.list, function(x) { x$Type }) == i) , ,drop=FALSE]

    if (nrow(category.GST.scores) > 60)
    {
      category.GST.order <- apply(category.GST.scores, 2, order, decreasing=TRUE)
      category.GST.order <- as.vector(t(category.GST.order))

      j <- 1

      while (length(unique(category.GST.order[1:(ncol(env$indata)*j)])) < 60)
      {
        j <- j + 1
      }

      top.scores <- category.GST.scores[unique(category.GST.order[1:(ncol(env$indata)*j)]),]
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
      paste(env$LETTERS[apply(summary.spots.fisher.p[rownames(top.scores), ,drop=FALSE], 1, which.min)],
            rownames(top.scores), sep="  ")

    colnames(top.scores) <- colnames(env$indata)

    heatmap(x=top.scores, cex.main=2,
                 col=env$color.palette.heatmaps(1000),
                 margins=c(10,20), scale="n", zlim=max(max(top.scores),-min(top.scores))*c(-1,1),
                 ColSideColors=env$group.colors, cexDend=0.6)

    par(new=TRUE, mar=c(3.5,29,35.5,2))

    image(matrix(c(1:1000), 1000, 1), axes=FALSE,
          col=env$color.palette.heatmaps(1000))

    box()
    axis(1, c(0,0.5,1), round(max(max(top.scores),-min(top.scores))*c(-1,0,1)), cex.axis=1.4)
    mtext("GSZ score", cex=2, line=32)
    mtext(paste("Category",i), cex=1, line=30)

    heatmap(x=top.scores, cex.main=2,
                 col=env$color.palette.heatmaps(1000),
                 margins=c(10,20), scale="n",  zlim=max(max(top.scores),-min(top.scores))*c(-1,1),
                 ColSideColors=env$group.colors, cexDend=0.6, Colv=NA)

    par(new=TRUE, mar=c(3.5,29,35.5,2))

    image(matrix(c(1:1000), 1000, 1), axes=FALSE,
          col=env$color.palette.heatmaps(1000))

    box()
    axis(1, c(0,0.5,1), round(max(max(top.scores),-min(top.scores))*c(-1,0,1)), cex.axis=1.4)
    mtext("GSZ score", cex=2, line=33.5)
    mtext(paste("Category",i), cex=1, line=32)

    o <- unlist(sapply(unique(env$group.labels), function(gr)
    {
      idx <- names(env$group.labels)[which(env$group.labels == gr)]

      if (length(idx) > 1)
      {
        hc <- hclust(dist(t(top.scores[,idx])))
        return(hc$labels[hc$order])
      }
      return(idx)
    }))

    heatmap(x=top.scores[,o], cex.main=2,
                 col=env$color.palette.heatmaps(1000),
                 margins=c(10,20), scale="n",
                 zlim=max(max(top.scores),-min(top.scores))*c(-1,1),
                 ColSideColors=env$group.colors[o], cexDend=0.6, Colv=NA)

    par(new=TRUE, mar=c(3.5,29,35.5,2) )

    image(matrix(c(1:1000), 1000, 1), axes=FALSE,
          col=env$color.palette.heatmaps(1000))

    box()
    axis( 1, c(0,0.5,1), round(max(max(top.scores),-min(top.scores))*c(-1,0,1)), cex.axis=1.4 )
    mtext( "GSZ score", cex=2, line=33.5 )
    mtext(paste("Category",i), cex=1, line=32)
  }

  dev.off()
}
