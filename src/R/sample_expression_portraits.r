pipeline.sampleExpressionPortraits <- function()
{
  ## Expression portraits
  filename <- file.path(paste(files.name, "- Results"), "Expression Portraits.pdf")
  util.info("Writing:", filename)

  pdf(filename, 29.7/2.54, 21/2.54)
  par(mfrow=c(7, 12))
  par(mar=c(0.3, 0.9, 4.5, 0.9))
  count.col <- 0

  for (gl in seq_along(unique(group.labels)))
  {
    plot(0, type="n", axes=FALSE, xlab="", ylab="", xlim=c(0,1))

    if (length(unique(group.labels)) > 1)
    {
      mtext(unique(group.labels)[gl], side=3, line=2, cex=1.5, at=0, font=3,
            adj=0, col=groupwise.group.colors[gl])
    }

    par(new=TRUE)

    for (j in which(group.labels == unique(group.labels)[gl]))
    {
      image(matrix(metadata[,j], preferences$dim.1stLvlSom, preferences$dim.1stLvlSom),
            axes=FALSE, col=colramp(1000))

      title(paste(j,":",colnames(indata)[j]), line=1, cex.main=0.8)

      if (length(unique(group.labels)) > 1)
      {
        title(paste("S =",round(group.silhouette.coef[j],2)), line=0.2,
              cex.main=0.6, col.main=groupwise.group.colors[gl])
      }

      box()
      count.col <- count.col + 1
    }

    if (count.col %% 12 != 0)
    {
      for (j in 1:(12 - count.col %% 12))
      {
        plot(0, type="n", axes=FALSE, xlab="", ylab="")
        count.col <- count.col + 1
      }
    }
  }

  par(mfrow=c(7, 12))
  par(mar=c(0.3, 0.9, 4.5, 0.9))
  count.col <- 0

  for (gl in seq_along(unique(group.labels)))
  {
    plot(0, type="n", axes=FALSE, xlab="", ylab="", xlim=c(0,1))

    if (length(unique(group.labels)) > 1)
    {
      mtext(unique(group.labels)[gl], side=3, line=2, cex=1.5, at=0, font=3,
            adj=0, col=groupwise.group.colors[gl])
    }

    par(new=TRUE)

    for (j in which(group.labels == unique(group.labels)[gl]))
    {
      md <- matrix(metadata[,j], preferences$dim.1stLvlSom, preferences$dim.1stLvlSom)
      nrz <- nrow(md)
      zfacet <- md[-1, -1] + md[-1, -nrz] + md[-nrz, -1] + md[-nrz, -nrz]
      facetcol <- cut(zfacet, 1000)

      persp(md, axes=FALSE, border=NA, expand=0.5, col=colramp(1000)[facetcol],
            phi=45, theta=-5, xlab="", ylab="", zlab="", box=TRUE)

      title(paste(j, ":", colnames(indata)[j]), line=1, cex.main=0.8)
      count.col <- count.col + 1
    }

    if (count.col %% 12 != 0)
    {
      for (j in 1:(12 - count.col %% 12))
      {
        plot(0, type="n", axes=FALSE, xlab="", ylab="")
        count.col <- count.col + 1
      }
    }
  }

  dev.off()


  ## Alternative expression portraits

  filename <- file.path(paste(files.name, "- Results"), "Expression Portraits - alternative scales.pdf")
  util.info("Writing:", filename)
  
  pdf(filename, 29.7/2.54, 21/2.54)
  
  
  # Absolute Metagene Portraits
  par(mar=c(0, 0, 0, 0), mfrow=c(1, 1))
  plot(0, type="n", xlab="", ylab="", axes=TRUE)
  text(1, 0.3, "Sample portraits in absolute color scale for all samples", cex=1.6)

  par(mfrow=c(7, 12))
  par(mar=c(0.3, 0.9, 4.5, 0.9))
  count.col <- 0

  for (gl in seq_along(unique(group.labels)))
  {
    plot(0, type="n", axes=FALSE, xlab="", ylab="", xlim=c(0, 1))

    if (length(unique(group.labels)) > 1)
    {
      mtext(unique(group.labels)[gl], side=3, line=2, cex=1.5, at=0, font=3,
            adj=0, col=groupwise.group.colors[gl])
    }

    par(new=TRUE)

    for (j in which(group.labels == unique(group.labels)[gl]))
    {
      image(matrix(metadata[,j], preferences$dim.1stLvlSom, preferences$dim.1stLvlSom),
            axes=FALSE, col = colramp(1000), zlim=c(min(metadata),max(metadata)))

      title(paste(j,":", colnames(indata)[j]), line=1, cex.main=0.8)
      box()
      count.col <- count.col + 1
    }

    if (count.col %% 12 != 0)
    {
      for (j in 1:(12 - count.col %% 12))
      {
        plot(0, type="n", axes=FALSE, xlab="", ylab="")
        count.col <- count.col + 1
      }
    }
  }


  # WAD Portraits
  par(mar=c(0, 0, 0, 0), mfrow=c(1, 1))
  plot(0, type="n", xlab="", ylab="", axes=TRUE)
  text(1, 0.3, "Sample portraits in WAD color scale", cex=1.6)

  par(mfrow=c(7, 12))
  par(mar=c(0.3, 0.9, 4.5, 0.9))
  count.col <- 0

  for (gl in seq_along(unique(group.labels)))
  {
    plot(0, type="n", axes=FALSE, xlab="", ylab="", xlim=c(0,1))

    if (length(unique(group.labels)) > 1)
    {
      mtext(unique(group.labels)[gl], side=3, line=2, cex=1.5, at=0, font=3,
            adj=0, col=groupwise.group.colors[gl])
    }

    par(new=TRUE)

    for (j in which(group.labels == unique(group.labels)[gl]))
    {
      image(matrix(WAD.metadata[,j], preferences$dim.1stLvlSom, preferences$dim.1stLvlSom),
            axes=FALSE, col=colramp(1000), cex.main=0.6)

      title(paste(j,":",colnames(indata)[j]), line=1, cex.main=0.8)
      box()
      count.col <- count.col + 1
    }

    if (count.col %% 12 != 0)
    {
      for (j in 1:(12 - count.col %% 12))
      {
        plot(0, type="n", axes=FALSE, xlab="", ylab="")
        count.col <- count.col + 1
      }
    }
  }

  # loglog FC Portraits
  par(mar=c(0, 0, 0, 0), mfrow=c(1, 1))
  plot(0, type="n", xlab="", ylab="", axes=TRUE)
  text(1, 0.3, "Sample portraits in loglog foldchange color scale", cex=1.6)

  par(mfrow=c(7, 12))
  par(mar=c(0.3, 0.9, 4.5, 0.9))
  count.col <- 0

  for (gl in seq_along(unique(group.labels)))
  {
    plot(0, type="n", axes=FALSE, xlab="", ylab="", xlim=c(0, 1))

    if (length(unique(group.labels)) > 1)
    {
      mtext(unique(group.labels)[gl], side=3, line=2, cex=1.5, at=0, font=3,
            adj=0, col=groupwise.group.colors[gl])
    }

    par(new=TRUE)

    for (j in which(group.labels == unique(group.labels)[gl]))
    {
      image(matrix(loglog.metadata[,j], preferences$dim.1stLvlSom, preferences$dim.1stLvlSom),
            axes=FALSE, col = colramp(1000))

      title(paste(j,":",colnames(indata)[j]), line=1, cex.main=0.8)
      box()
      count.col <- count.col + 1
    }

    if (count.col %% 12 != 0)
    {
      for (j in 1:(12 - count.col %% 12))
      {
        plot(0, type="n", axes=FALSE, xlab="", ylab="")
        count.col <- count.col + 1
      }
    }
  }

  dev.off()
}
