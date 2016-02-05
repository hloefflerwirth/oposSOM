sort.label <- function(x)
{
  numeric <- function(x)
  {
    suppressWarnings(as.numeric(x))
  }

  nonnumeric <- function(x)
  {
    ifelse(is.na(suppressWarnings(as.numeric(x))), toupper(x), NA)
  }

  x <- as.character(x)
  delim <- "\\$\\@\\$"
  delimited <- gsub("([+-]{0,1}[0-9]+([eE][\\+\\-]{0,1}[0-9]+){0,1})",
                    paste(delim,"\\1",delim,sep=""), x)

  step1 <- strsplit(delimited, delim)
  step1 <- lapply(step1, function(x) { x[x>""] })

  step1.numeric <- lapply(step1, numeric)
  step1.character <- lapply(step1, nonnumeric)

  maxelem <- max(sapply(step1, length))

  step1.numeric.t <- lapply(1:maxelem, function(i)
  {
    sapply(step1.numeric, function(x) { x[i] })
  })

  step1.character.t <- lapply(1:maxelem, function(i)
  {
    sapply(step1.character, function(x) { x[i] })
  })

  rank.numeric <- sapply(step1.numeric.t,rank)
  rank.character <- sapply(step1.character.t, function(x) { as.numeric(factor(x)) })
  rank.numeric[!is.na(rank.character)] <- 0
  rank.character <- t(t(rank.character) + apply(matrix(rank.numeric),2,max,na.rm=TRUE))
  rank.overall <- ifelse(is.na(rank.character), rank.numeric,rank.character)

  order.frame <- as.data.frame(rank.overall)
  x <- x[do.call("order", order.frame)]
  return(x)
}

pipeline.chromosomeExpressionReports <- function()
{
  if(ncol(metadata) > 1000 || length(gene.positions.list) == 0 ) return()
  
  # prepare chromosome geneset lists
  chr.gs.list <- lapply(gene.positions.list, function(x)
  {
    list(Genes=gene.ids[unlist(x)], Type="Chr")
  })

  names(chr.gs.list) <- paste("Chr", names(gene.positions.list))
  chr.pq.gs.list <- list()

  for (i in seq_along(gene.positions.list))
  {
    genes <- unlist(gene.positions.list[[i]])
    names(genes) <- substr(names(genes),1,1)

    pq.list <- tapply(genes, names(genes), c)
    names(pq.list) <- paste("Chr", names(gene.positions.list)[i], names(pq.list))

    chr.pq.gs.list <- c(chr.pq.gs.list, pq.list)
  }

  chr.pq.gs.list <- lapply(chr.pq.gs.list, function(x)
  {
    list(Genes=gene.ids[x], Type="Chr")
  })

  
  # Calculate GSZ
  util.progress(0, 10)

  sample.chr.GSZ <- do.call(cbind, lapply( 1:ncol(indata), function(m)
  {
    return(GeneSet.GSZ(unique.protein.ids, t.ensID.m[,m], chr.gs.list, sort=FALSE))
  }))

  dimnames(sample.chr.GSZ) <- list(names(chr.gs.list), colnames(indata))

  util.progress(5, 10)

  sample.chr.pq.GSZ <- do.call(cbind, lapply( 1:ncol(indata), function(m)
  {
    return(GeneSet.GSZ(unique.protein.ids, t.ensID.m[,m], chr.pq.gs.list, sort=FALSE))
  }))

  dimnames(sample.chr.pq.GSZ) <- list(names(chr.pq.gs.list), colnames(indata))

  util.progress.terminate()


  
  # Heatmap Outputs
  filename <- file.path(paste(files.name, "- Results"),
                        "Geneset Analysis",
                        "0verview Chromosome Expression.pdf")

  util.info("Writing:", filename)
  pdf(filename, 21/2.54, 21/2.54)

  heatmap.wrap(x=sample.chr.GSZ, cex.main=2,
               col=color.palette.heatmaps(1000),
               mar=c(10,20), scale="n",
               zlim=max(max(sample.chr.GSZ),-min(sample.chr.GSZ))*c(-1,1),
               ColSideColors=group.colors, cexDend=0.6)

  par(new=TRUE, mar=c(3.5,29,35.5,2))

  image(matrix(c(1:1000), 1000, 1), axes=FALSE,
        col=color.palette.heatmaps(1000))

  box()

  axis(1, c(0,0.5,1),
       round(max(max(sample.chr.GSZ),-min(sample.chr.GSZ))*c(-1,0,1)), cex.axis=1)

  title(main="GSZ score", line=0.5, cex.main=0.8)

  heatmap.wrap(x=sample.chr.GSZ, cex.main=2,
               col=color.palette.heatmaps(1000),
               mar=c(10,20), scale="n", zlim=max(max(sample.chr.GSZ),-min(sample.chr.GSZ))*c(-1,1),
               ColSideColors=group.colors, cexDend=0.6, Colv=NA)

  par(new=TRUE, mar=c(3.5,29,35.5,2))

  image(matrix(c(1:1000), 1000, 1), axes=FALSE,
        col=color.palette.heatmaps(1000))

  box()

  axis(1, c(0,0.5,1),
       round(max(max(sample.chr.GSZ),-min(sample.chr.GSZ))*c(-1,0,1)), cex.axis=1)

  title(main="GSZ score", line=0.5, cex.main=0.8)

  heatmap.wrap(x=sample.chr.GSZ[rev(sort.label(rownames(sample.chr.GSZ))),], cex.main=2,
               col=color.palette.heatmaps(1000),
               mar=c(10,20), scale="n", zlim=max(max(sample.chr.GSZ),-min(sample.chr.GSZ))*c(-1,1),
               ColSideColors=group.colors, cexDend=0.6, Colv=NA, Rowv=NA)

  par(new=TRUE, mar=c(3.5,29,35.5,2))

  image(matrix(c(1:1000), 1000, 1), axes=FALSE,
        col=color.palette.heatmaps(1000))

  box()

  axis(1, c(0,0.5,1),
       round(max(max(sample.chr.GSZ),-min(sample.chr.GSZ))*c(-1,0,1)), cex.axis=1)

  title(main="GSZ score", line=0.5, cex.main=0.8)

  heatmap.wrap(x=sample.chr.pq.GSZ, cex.main=2,
               col=color.palette.heatmaps(1000),
               mar=c(10,20), scale="n", zlim=max(max(sample.chr.pq.GSZ),-min(sample.chr.pq.GSZ))*c(-1,1),
               ColSideColors=group.colors, cexDend=0.6)

  par(new=TRUE, mar=c(3.5,29,35.5,2))

  image(matrix(c(1:1000), 1000, 1), axes=FALSE,
        col=color.palette.heatmaps(1000))

  box()

  axis(1, c(0,0.5,1),
       round(max(max(sample.chr.pq.GSZ),-min(sample.chr.pq.GSZ))*c(-1,0,1)), cex.axis=1)

  title(main="GSZ score", line=0.5, cex.main=0.8)

  heatmap.wrap(x=sample.chr.pq.GSZ, cex.main=2,
               col=color.palette.heatmaps(1000),
               mar=c(10,20), scale="n", zlim=max(max(sample.chr.pq.GSZ),-min(sample.chr.pq.GSZ))*c(-1,1),
               ColSideColors=group.colors, cexDend=0.6, Colv=NA)

  par(new=TRUE, mar=c(3.5,29,35.5,2))

  image(matrix(c(1:1000), 1000, 1), axes=FALSE,
        col=color.palette.heatmaps(1000))

  box()

  axis(1, c(0,0.5,1),
       round(max(max(sample.chr.pq.GSZ),-min(sample.chr.pq.GSZ))*c(-1,0,1)), cex.axis=1)

  title(main="GSZ score", line=0.5, cex.main=0.8)


  heatmap.wrap(x=sample.chr.pq.GSZ[rev(sort.label(rownames(sample.chr.pq.GSZ))),],
               cex.main=2, col=color.palette.heatmaps(1000),
               mar=c(10,20), scale="n", zlim=max(max(sample.chr.pq.GSZ),-min(sample.chr.pq.GSZ))*c(-1,1),
               ColSideColors=group.colors, cexDend=0.6, Colv=NA, Rowv=NA)

  par(new=TRUE, mar=c(3.5,29,35.5,2))

  image(matrix(c(1:1000), 1000, 1), axes=FALSE,
        col=color.palette.heatmaps(1000))

  box()
  axis(1, c(0,0.5,1),
       round(max(max(sample.chr.pq.GSZ),-min(sample.chr.pq.GSZ))*c(-1,0,1)), cex.axis=1)

  title(main="GSZ score", line=0.5, cex.main=0.8)

  dev.off()
}
