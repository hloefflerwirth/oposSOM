## Heatmap Wrapper cause bad Parameter-Handling Row/ColSideColors
heatmap.wrap <- function(RowSideColors, ColSideColors, ...)
{
  if (!missing(RowSideColors) && !missing(ColSideColors))
  {
    if (length(unique(RowSideColors)) > 1 && length(unique(ColSideColors)) > 1)
    {
      heatmap(RowSideColors=RowSideColors, ColSideColors=ColSideColors, ...)
      return(NULL)
    }
  } else if (!missing(ColSideColors))
  {
    if (length(unique(ColSideColors)) > 1)
    {
      heatmap(ColSideColors=ColSideColors, ...)
      return(NULL)
    }
  }

  heatmap(...)
}

pipeline.2ndLvlSimilarityAnalysis <- function()
{

  # Filter Metagenes
  metagene.filter.list <<- list()

  metagene.filter.list[[1]] <<- list(s=c(1:preferences$dim.1stLvlSom^2),
                                     n=paste("all", preferences$dim.1stLvlSom^2, "Metagenes"))

  spot.metagenes <- unique(unlist(sapply(spot.list.overexpression$spots, function(x) { x$metagenes })))

  metagene.filter.list[[2]] <<- list(s=spot.metagenes,
                                     n=paste(length(spot.metagenes), "Spot-Metagenes"))

  if (ncol(metadata) < 1000 && preferences$dim.1stLvlSom^2 > 100)
  {
    metagene.filter.list[[length(metagene.filter.list)+1]] <<-
      list(s=order(apply(abs(metadata), 1, max), decreasing=TRUE)[1:100],
           n="High Expression: 100 Metagenes")
  }

  if (ncol(metadata) < 1000 && preferences$dim.1stLvlSom^2 > 100)
  {
    metagene.filter.list[[length(metagene.filter.list)+1]] <<-
      list(s=order(apply(abs(metadata), 1, var), decreasing=TRUE)[1:100],
           n="Variance: 100 Metagenes")
  }

  filename <- file.path(paste(files.name, "- Results"), "2nd lvl Metagene Analysis", "Similarity Analysis.pdf")
  util.info("Writing:", filename)
  pdf(filename, 29.7/2.54, 21/2.54)

  for (i in seq_along(metagene.filter.list))
  {
    s <- metadata[metagene.filter.list[[i]]$s ,]
    par(mar=c(1,1,1,1))

    if (ncol(metadata) > 2)
    {
      phylo.tree <- nj(dist(t(s)))
      phylo.tree$tip.label <- colnames(indata)
      plot.phylo(phylo.tree, "unrooted", cex=0.5, tip.color=group.colors)
      title(paste("Phylogenetic sample clustering,",metagene.filter.list[[i]]$n))
      box()
      
      par(new=TRUE)
      plot(0,type="n", axes=FALSE, xlab="", ylab="")    
      legend("bottomright", as.character(unique(group.labels)), cex=0.5,
             text.col=groupwise.group.colors, bg="white")
    }

    heatmap.wrap(x=s, col=colramp(1000), main=paste("Clustering heatmap,",metagene.filter.list[[i]]$n),
                 margins=c(10, 5), scale="n", labRow=NA, ColSideColors=group.colors)

    par(new=TRUE)
    plot(0,type="n", axes=FALSE, xlab="", ylab="")

    legend("bottomright", as.character(unique(group.labels)), cex=0.5,
           text.col=groupwise.group.colors, bg="white")
    
    par(new=TRUE, mar = c(25, 55, 10.8, 2))
    image(matrix(1:100, 1, 100), col = colramp(1000), axes=FALSE)
    axis(2, round(range(s),1), at=c(0, 1), las=2, tick=FALSE, pos=0, cex.axis=1)
    

    heatmap.wrap(x=s, col=colramp(1000), main=paste("Clustering heatmap,",metagene.filter.list[[i]]$n),
                 margins=c(10, 6), scale="n", labRow=NA, ColSideColors=group.colors, Colv=NA)

    par(new=TRUE)
    plot(0,type="n", axes=FALSE, xlab="", ylab="")

    legend("bottomright", as.character(unique(group.labels)), cex=0.5,
           text.col=groupwise.group.colors, bg="white")
    
    par(new=TRUE, mar = c(31.6, 55, 4.2, 2))
    image(matrix(1:100, 1, 100), col = colramp(1000), axes=FALSE)
    axis(2, round(range(s),1), at=c(0, 1), las=2, tick=FALSE, pos=0, cex.axis=1)
    
  }

  dev.off()
}
