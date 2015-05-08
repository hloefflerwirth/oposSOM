pipeline.2ndLvlComponentAnalysis <- function()
{
  filename <- file.path(paste(files.name, "- Results"), "2nd lvl Metagene Analysis", "Component Analysis.pdf")
  util.info("Writing:", filename)

  pdf(filename, 21/2.54, 29.7/2.54, useDingbats=FALSE)

  for (i in seq_along(metagene.filter.list))
  {
    try({
      ICA.metagenes <- fastICA(t(metadata[metagene.filter.list[[i]]$s ,]), 3)$S
      z <- ICA.metagenes[,2]

      layout(matrix(c(0, 1, 0)),heights=c(1, 3, 1))
      par(mar=c(1, 1, 1, 1))

      scatterplot3d(ICA.metagenes,
                    cex.symbols=4*(1-((z-min(z))/(max(z)-min(z))))+2,
                    color=group.colors,
                    pch=16, tick.marks=FALSE, xlab="", ylab="", zlab="",
                    main=paste("Independent component analysis,",metagene.filter.list[[i]]$n), mar=c(1,1,1,1))

      par(new=TRUE)
      plot(0,type="n",axes=FALSE,xlab="",ylab="")
      text(0.75, 0.7485, "Component 2", cex=1, srt=38)
      mtext("component 1", 1, cex=0.8, line=-1, at=0.84)
      mtext("component 3", 2, cex=0.8, line=-1, at=-0.3)
      par(new=TRUE)

      scatterplot3d(ICA.metagenes,
                    cex.symbols=4*(1-((z-min(z))/(max(z)-min(z))))+2,
                    color="black", pch=1, tick.marks=FALSE, xlab="", ylab="",
                    zlab="", axis=FALSE, grid=FALSE, mar=c(1,1,1,1))

      par(new=TRUE)
      plot(0,type="n", axes=FALSE, xlab="", ylab="")

      legend("bottomright", as.character(unique(group.labels)), cex=0.5,
             text.col=groupwise.group.colors, bg="white")

      
      
      layout(matrix(c(1,2)))
      par(mar=c(0.1,3,1,3))

      plot(ICA.metagenes[,1], ICA.metagenes[,3], type="p", pch=16,
           col=group.colors, cex=3, axes=FALSE, xlab="", ylab="",
           main=paste("Independent component analysis,",metagene.filter.list[[i]]$n), cex.main=0.8)

      mtext("component 3",2,cex=0.8)
      points(ICA.metagenes[,1], ICA.metagenes[,3], pch=16, col=group.colors, cex=3)
      points(ICA.metagenes[,1], ICA.metagenes[,3], pch=1, col="black", cex=3)
      text(ICA.metagenes[,1], ICA.metagenes[,3], colnames(indata), col="gray20", cex=0.6)
      box()

      par(mar=c(1,3,0.1,3))

      plot(ICA.metagenes[,1], ICA.metagenes[,2], type="p", pch=16,
           col=group.colors, cex=3, axes=FALSE, xlab="", ylab="", main="")

      mtext("component 1",1,cex=0.8)
      mtext("component 2",2,cex=0.8)
      
      points(ICA.metagenes[,1], ICA.metagenes[,2], pch=16, col=group.colors, cex=3)
      points(ICA.metagenes[,1], ICA.metagenes[,2], pch=1, col="black", cex=3)
      text(ICA.metagenes[,1], ICA.metagenes[,2], colnames(indata), col="gray20", cex=0.6)
      box()
      
      if (length(unique(group.labels)) > 1)
      {
        legend("topright", as.character(unique(group.labels)), cex=0.4,
               text.col=groupwise.group.colors, bg="white")
      }

    }, silent=TRUE)
  }

  dev.off()
}
