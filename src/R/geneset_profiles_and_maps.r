pipeline.genesetProfilesAndMaps <- function()
{
  samples.GSZ.scores <<- do.call(cbind, lapply(spot.list.samples, function(x)
  {
    x$GSZ.score[names(gs.def.list)]
  }))

  ### CSV output ###
  filename <- file.path(paste(files.name, "- Results"), "CSV Sheets", "Sample GSZ scores.csv")
  util.info("Writing:", filename)
  write.csv2(samples.GSZ.scores, filename)

  progress.current <- 0
  progress.max <- nrow(samples.GSZ.scores) + length(gs.def.list)
  util.progress(progress.current, progress.max)

  ## Geneset Profiles over Samples
  if (preferences$geneset.analysis.exact)
  {
    fdr.threshold <- 0.1

    GSZs <-
      cbind(as.vector(unlist(sapply(spot.list.samples, function(x) { x$GSZ.score[names(gs.def.list)] }))),
            as.vector(unlist(sapply(spot.list.samples, function(x) { x$GSZ.p.value[names(gs.def.list)] }))))

    GSZs <- GSZs[which(!is.na(GSZs[,2])),]

    fdrtool.result <- fdrtool(GSZs[,2], statistic="pvalue", verbose=FALSE, plot=FALSE)
    fdr.significant.spot.list.samples <- length(which(fdrtool.result$lfdr < fdr.threshold))
    fdr.gsz.threshold <- sort(GSZs[,1], decreasing=TRUE)[fdr.significant.spot.list.samples]
  } else
  {
    fdr.gsz.threshold <- 0
  }

  dirname <- file.path(paste(files.name, "- Results"), "Geneset Analysis")
  util.info("Writing:", file.path(dirname, "*.{csv,pdf}"))

  for (i in 1:nrow(samples.GSZ.scores))
  {
    filename.prefix <- substring(make.names(names(gs.def.list)[i]), 1, 100)
    pdf(file.path(dirname, paste(filename.prefix, "profile.pdf")), 29.7/2.54, 21/2.54)

    par(mar=c(15,6,4,5))

    leading.metagenes <-
      as.numeric(names(sort(table(som.nodes[names(gene.ids)[which(gene.ids %in% gs.def.list[[i]]$Genes)]]),
                            decreasing=TRUE))[1:5])

    ylim <- c(-10, 20)

    bar.coords <- barplot(samples.GSZ.scores[i,], beside=TRUE,
                          main=rownames(samples.GSZ.scores)[i], las=2,
                          cex.names=1.2, col=group.colors, cex.main=1,
                          ylim=ylim,
                          border=if (ncol(indata) < 80) "black" else NA ,
                          names.arg=if (ncol(indata)<100) colnames(indata) else rep("",ncol(indata)))

    abline(h=fdr.gsz.threshold, lty=2)
    abline(h=-fdr.gsz.threshold, lty=2)

    mtext("GSZ", side=2, line=3.5, cex=2)

    mean.boxes <- by(samples.GSZ.scores[i,], group.labels, c)[unique(group.labels)]

    boxplot(mean.boxes, col=groupwise.group.colors, main=rownames(samples.GSZ.scores)[i],
            cex.main=1, ylim=ylim, axes=FALSE, yaxs="i")

    axis(1, seq_along(unique(group.labels)), unique(group.labels), las=2, tick=FALSE)
    axis(2, las=2)

    abline(h=0, lty=2)
    abline(h=fdr.gsz.threshold, lty=2)
    abline(h=-fdr.gsz.threshold, lty=2)

    mtext("GSZ", side=2, line=3.5, cex=2)

    #################################################

    spot.fisher.p <- -log10(sapply(spot.list.overexpression$spots, function(x)
    {
      x$Fisher.p[names(gs.def.list)[i]]
    }))

    names(spot.fisher.p) <- LETTERS[seq_along(spot.fisher.p)]

    barplot(spot.fisher.p, main=paste("Enrichment in spots:", names(gs.def.list)[i]),
            las=1, cex.names=1.2)

    mtext("log( p.value )", side=2, line=3.5, cex=2)

    #################################################

    spot.expression <- t(sapply(spot.list.overexpression$spots, function(x)
    {
      g <- intersect(gene.ids[x$genes], gs.def.list[[i]]$Genes)
      g <- names(gene.ids)[which(gene.ids %in% g)]

      if (length(g) > 0)
      {
        return(apply(indata[g,,drop=FALSE], 2, mean))
      } else
      {
        return(rep(0, ncol(indata)))
      }
    }))

    colnames(spot.expression) <- colnames(indata)

    for (ii in 1:nrow(spot.expression))
    {
      barplot(spot.expression[ii,], beside=TRUE,
              main=paste("Expression of",names(gs.def.list)[i],"in Spot",LETTERS[ii]),
              las=2, cex.names=1.2, col=group.colors, cex.main=1, ylim=range(spot.expression),
              border=if (ncol(indata) < 80) "black" else NA)

      mtext(expression(paste(Delta,"e")), side=2, line=3.5, cex=2)
    }

    dev.off()

    #################################################

    genes <- names(gene.ids)[which(gene.ids %in% gs.def.list[[i]]$Genes)]

    out <- data.frame(AffyID=names(gene.ids[genes]),
                      EnsemblID=gene.ids[genes],
                      Metagene=gene.coordinates[genes],
                      Max.expression.sample=colnames(indata)[apply(indata[genes, ,drop=FALSE], 1, which.max)],
                      GeneSymbol=gene.names[genes],
                      Description=gene.descriptions[genes])

    write.csv2(out, file.path(dirname, paste(filename.prefix, ".csv", sep="")), row.names=FALSE)

    progress.current <- progress.current + 1
    util.progress(progress.current, progress.max)
  }

  ## Gensets Population Maps
  for (i in seq_along(gs.def.list))
  {
    filename <- paste(substring(make.names(names(gs.def.list)[i]), 1, 100), "map.pdf")
    pdf(file.path(dirname, filename), 21/2.54, 21/2.54)

    ### Population Map ###
    n.map <- matrix(0,preferences$dim.1stLvlSom,preferences$dim.1stLvlSom)
    gs.nodes <- som.nodes[names(gene.ids)[which(gene.ids %in% gs.def.list[[i]]$Genes)]]
    n.map[as.numeric(names(table(gs.nodes)))] <- table(gs.nodes)
    n.map[which(n.map==0)] <- NA
    n.map <- matrix(n.map, preferences$dim.1stLvlSom)

    x.test <-
      sum((apply(n.map, 1, sum, na.rm=TRUE) - sum(n.map, na.rm=TRUE)*(1/preferences$dim.1stLvlSom))^2 /
          sum(n.map, na.rm=TRUE) *
          (1/preferences$dim.1stLvlSom))

    y.test <-
      sum((apply(n.map, 2, sum, na.rm=TRUE) - sum(n.map, na.rm=TRUE)*(1/preferences$dim.1stLvlSom))^2 /
          sum(n.map, na.rm=TRUE) *
          (1/preferences$dim.1stLvlSom))

    chi.sq.p.value <- 1-pchisq(x.test+y.test, 1)

    par(mfrow=c(1, 1))
    par(mar=c(5, 6, 4, 5))

    lim <- c(1,preferences$dim.1stLvlSom) + preferences$dim.1stLvlSom * 0.01 * c(-1, 1)

    colr <- colramp(1000)[(na.omit(as.vector(n.map)) - min(n.map,na.rm=TRUE)) /
                          max(1, (max(n.map,na.rm=TRUE) - min(n.map,na.rm=TRUE))) *
                          999 + 1]

    plot(which(!is.na(n.map), arr.ind=TRUE), xlim=lim, ylim=lim, pch=16, axes=FALSE,
         xlab="",ylab="", xaxs="i", yaxs="i", main=names(gs.def.list)[i], col=colr,
         cex=0.5 + na.omit(as.vector(n.map)) / max(n.map,na.rm=TRUE) * 2.8)

    title(sub=paste("# features =", sum(gene.ids %in% gs.def.list[[i]]$Genes)),line=0)
    title(sub=paste("chi-square p = ",round(chi.sq.p.value,2)),line=1)
    box()

    par(new=TRUE, mar=c(1,0,0,0))
    layout(matrix(c(0, 0, 0, 0, 1, 0, 0, 0, 0), 3, 3), c(1, 0.05, 0.02), c(0.15, 0.3, 1))
    image(matrix(1:100, 1, 100), col = colramp(1000), axes=FALSE)
    axis(2, at=c(0,1), c(min(n.map,na.rm=TRUE), max(n.map,na.rm=TRUE)), las=2, tick=FALSE, pos=-0.5)
    box()

    ### Smooth Population Map ###
    if (any(!is.na(n.map)))
    {
      par(mfrow=c(1, 1))
      par(mar=c(5, 6, 4, 5))

      coords <-
        som.result$visual[which(rownames(indata) %in% names(gene.ids[which(gene.ids %in% gs.def.list[[i]]$Genes)])),
                          c(1,2), drop=FALSE]

      if (nrow(coords) == 1)
      {
        coords <- rbind(coords, coords)
      }

      suppressWarnings(smoothScatter(coords+1 ,
                                     main=names(gs.def.list)[i],
                                     cex.main=2,
                                     xlim=c(1,preferences$dim.1stLvlSom),
                                     ylim=c(1,preferences$dim.1stLvlSom),
                                     xaxs="i", yaxs="i", axes=FALSE, xlab="",
                                     ylab="", nrpoints=0))

      title(sub=paste("# features =", sum(gene.ids %in% gs.def.list[[i]]$Genes),
                      ", max =", max(n.map,na.rm=TRUE)),line=0)
    }

    dev.off()

    progress.current <- progress.current + 1
    util.progress(progress.current, progress.max)
  }

  util.progress.terminate()
}
