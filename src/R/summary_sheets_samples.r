pipeline.summarySheetsSamples <- function()
{
  if(ncol(metadata) > 1000) return()
  
  dir.create(output.paths["Summary Sheets Samples"], showWarnings=FALSE)

  #### Summary Sheets ####
  n.genes.in.genesets <- 0

  if (preferences$geneset.analysis)
  {
    n.genes.in.genesets <- length(intersect(unique(unlist(gs.def.list)), gene.ids))
  }

  ylim.max <- 0

  for (m in 1:ncol(indata))
  {
    h <- hist(p.g.m[,m], bre=20, plot=FALSE)
    y.max <- max(h$density)

    if (y.max > ylim.max)
    {
      ylim.max <- y.max
    }
  }

  util.info("Writing:", file.path(output.paths["Summary Sheets Samples"], "*.pdf"))

  for (m in 1:ncol(indata))
  {
    basename <- paste(make.names(make.unique(colnames(indata))[m]), ".pdf", sep="")
    pdf(file.path(output.paths["Summary Sheets Samples"], basename), 29.7/2.54, 21/2.54)

    ## Global Sheet
    layout(matrix(c(1,2,4,1,3,0,5,5,6,7,7,8), 3, 4), widths=c(1,1,2,2), heights=c(2,1,1))
    par(mar=c(0,0,0,0))
    plot(0, type="n", axes=FALSE, xlab="", ylab="", xlim=c(0,1), ylim=c(0,1))

    text(0.1, 0.94, colnames(indata)[m] , cex=3, adj=0)
    text(0.1, 0.8, "Global Summary" , cex=1.8, adj=0)
    text(0.1, 0.7,  paste("%DE =", round(perc.DE.m[colnames(indata)[m]], 2)), adj=0)

    all.fdr.genes <- which(fdr.g.m[,m] < 0.2)
    plus.fdr.genes <- which(indata[all.fdr.genes, m] > 0)
    minus.fdr.genes <- which(indata[all.fdr.genes, m] <= 0)

    text(0.1, 0.65, paste("# genes with fdr < 0.2  =",
                          length(all.fdr.genes), " (",
                          length(plus.fdr.genes), "+ /",
                          length(minus.fdr.genes), " -)"), adj=0)

    all.fdr.genes <- which(fdr.g.m[,m] < 0.1)
    plus.fdr.genes <- which(indata[all.fdr.genes, m] > 0)
    minus.fdr.genes <- which(indata[all.fdr.genes, m] <= 0)

    text(0.1, 0.6, paste("# genes with fdr < 0.1  =",
                         length(all.fdr.genes)," (",
                         length(plus.fdr.genes), "+ /",
                         length(minus.fdr.genes), " -)"), adj=0)

    all.fdr.genes <- which(fdr.g.m[,m] < 0.05)
    plus.fdr.genes <- which(indata[all.fdr.genes, m] > 0)
    minus.fdr.genes <- which(indata[all.fdr.genes, m] <= 0)

    text(0.1, 0.55, paste("# genes with fdr < 0.05  =",
                          length(all.fdr.genes)," (",
                          length(plus.fdr.genes), "+ /",
                          length(minus.fdr.genes), " -)"), adj=0)

    all.fdr.genes <- which(fdr.g.m[,m] < 0.01)
    plus.fdr.genes <- which(indata[all.fdr.genes, m] > 0)
    minus.fdr.genes <- which(indata[all.fdr.genes, m] <= 0)

    text(0.1, 0.5, paste("# genes with fdr < 0.01 =",
                         length(all.fdr.genes)," (",
                         length(plus.fdr.genes), "+ /",
                         length(minus.fdr.genes), " -)"), adj=0)

    text(0.1, 0.425, paste("# genes in genesets =", n.genes.in.genesets), adj=0)
    text(0.1, 0.35,  paste("<FC> =", round(mean(indata[,m]), 2)), adj=0)
    text(0.1, 0.3, paste("<t-score> =", round(mean(t.g.m[,m]), 2)), adj=0)
    text(0.1, 0.25,  paste("<p-value> =", round(10 ^ mean(log10(p.g.m[,m])), 2)), adj=0)
    text(0.1, 0.2, paste("<fdr> =", round(mean(fdr.g.m[,m]), 2)), adj=0)

    par(mar=c(2,3,3,1))

    image(matrix(metadata[,m], preferences$dim.1stLvlSom, preferences$dim.1stLvlSom),
          axes=FALSE, col = colramp(1000), main="Profile", cex.main=1.5)

    axis(1,
         seq(0, 1, length.out=preferences$dim.1stLvlSom/10+1),
         c(1, seq(10, preferences$dim.1stLvlSom, length.out=preferences$dim.1stLvlSom/10)),
         cex.axis=1.0)

    axis(2,
         seq(0, 1, length.out=preferences$dim.1stLvlSom/10+1),
         c(1, seq(10, preferences$dim.1stLvlSom, length.out=preferences$dim.1stLvlSom/10)),
         cex.axis=1.0, las=1)

    box()

    image(matrix(metadata[,m], preferences$dim.1stLvlSom, preferences$dim.1stLvlSom),
          axes=FALSE, col=colramp(1000), main="Regulated Spots", cex.main=1.5)

    par(new=TRUE)

    mask <- spot.list.samples[[m]]$regulated
    mask[which(is.na(spot.list.samples[[m]]$regulated))] <- 1
    mask[which(!is.na(spot.list.samples[[m]]$regulated))] <- NA

    image(matrix(mask, preferences$dim.1stLvlSom, preferences$dim.1stLvlSom),
          axes=FALSE, col = "white")

    axis(1,
         seq(0, 1, length.out=preferences$dim.1stLvlSom/10+1),
         c(1, seq(10, preferences$dim.1stLvlSom, length.out=preferences$dim.1stLvlSom/10)),
         cex.axis=1.0)

    axis(2,
         seq(0, 1, length.out=preferences$dim.1stLvlSom/10+1),
         c(1, seq(10, preferences$dim.1stLvlSom, length.out=preferences$dim.1stLvlSom/10)),
         cex.axis=1.0, las=1)

    box()
    par(mar=c(2,8,3,6))
    plot(0, type="n", axes=FALSE, xlab="", ylab="", xlim=c(0,1), ylim=c(0,1))
    par(mar=c(0,0,0,0))

    n.genes <- 20
    x.coords <- c(0, 0.06, 0.2, 0.28, 0.36, 0.44, 0.52)
    y.coords <- seq(0.75, 0.02, length.out=n.genes)

    o <- order(p.g.m[,m])[1:n.genes]
    plot(0, type="n", axes=FALSE, xlab="", ylab="", xlim=c(0,1), ylim=c(0,1))

    text(0, 0.88, "Global Genelist", cex=1.8, adj=0)

    text(x.coords, rep(c(0.82, 0.80), 4)[1:7],
         c("Rank", "ID", "log(FC)", "p-value", "fdr", "Metagene", "Description"),
         cex=1, adj=0)

    text(x.coords[1], y.coords, c(1:n.genes), adj=0)
    text(x.coords[2], y.coords, rownames(indata)[o], cex=0.6, adj=0)
    rect(x.coords[3]-0.02, y.coords[1]+0.01, 1, 0, border="white", col="white")
    text(x.coords[3], y.coords, round(indata[o, m], 2), cex=0.6, adj=0)
    text(x.coords[4], y.coords, format(p.g.m[o, m], digits=1), cex=0.6, adj=0)
    text(x.coords[5], y.coords, format(fdr.g.m[o, m], digits=1), cex=0.6, adj=0)
    text(x.coords[6], y.coords, gene.coordinates[o], cex=0.6, adj=0)
    text(x.coords[7], y.coords, gene.descriptions[o], cex=0.6, adj=0)

    par(mar=c(3,6,2,6))

    hist(p.g.m[,m], bre=20, freq=FALSE, xlab="p-value", ylab="", main="p-values",
         ylim=c(0,ylim.max), las=1, cex.main=1.5, cex.lab=1, cex.axis=1)

    box()
    mtext("Density", side=2, line=3, cex=1)
    mtext("FDR", side=4, line=3, cex=1)
    mtext(paste("%DE =", round(perc.DE.m[colnames(indata)[m]] ,2)), line=-1.2, cex=0.5)

    abline(h=n.0.m[colnames(indata)[m]], col="gray", lwd=2)

    par(new=TRUE)
    plot(0, type="n", xlim=c(0,1), ylim=c(0,1), xlab="", ylab="", axes=FALSE)
    axis(4, seq(0, 1, 0.2), seq(0, 1, 0.2), las=1, cex.axis=1)
    o <- order(p.g.m[,m])
    lines(p.g.m[o,m], Fdr.g.m[o,m], lty=2, lwd=2)
    lines(p.g.m[o,m], fdr.g.m[o,m], lty=3, lwd=3)

    legend("topright", c("p", expression(eta[0]), "Fdr", "fdr"),
           col=c("black","gray","black","black"), lty=c(1,1,2,3),
           lwd=c(1,1,1,2), cex=0.7)

    if (preferences$geneset.analysis)
    {
      n.sets <- 20

      top.gs.score <- sort(spot.list.samples[[m]]$GSZ.score, decreasing=TRUE)[1:n.sets]
      top.gs.p <- spot.list.samples[[m]]$GSZ.p.value[names(top.gs.score)]

      par(mar=c(0,0,0,0))

      x.coords <- c(0, 0.1, 0.18, 0.30, 0.39, 0.47)
      y.coords <- seq(0.75, 0.4, length.out=n.sets)

      plot(0, type="n", axes=FALSE, xlab="", ylab="", xlim=c(0,1), ylim=c(0,1))

      text(0, 0.88, "Global Geneset Analysis", cex=1.8, adj=0)
      text(x.coords, 0.82, c("Rank", "GSZ", "p-value", "#all", "Geneset", ""), cex=1, adj=0)
      text(x.coords[1], 0.77, "Overexpressed", cex=0.8, adj=0, font=3)

      text(x.coords[1], y.coords, c(1:n.genes), adj=0)
      text(x.coords[2], y.coords, round(top.gs.score, 2), cex=0.6, adj=0)
      text(x.coords[3], y.coords, format(top.gs.p, digits=1), cex=0.6, adj=0)

      text(x.coords[4], y.coords, sapply(gs.def.list[names(top.gs.score)],
                                         function(x) { length(x$Genes) }), cex=0.6, adj=0)

      text(x.coords[5], y.coords, gs.def.list.categories[names(top.gs.score)], cex=0.6, adj=0)
      text(x.coords[6], y.coords, names(top.gs.score), cex=0.6, adj=0)

      top.gs.score <- sort(spot.list.samples[[m]]$GSZ.score, decreasing=FALSE)[1:n.sets]
      top.gs.p <- spot.list.samples[[m]]$GSZ.p.value[names(top.gs.score)]

      y.coords <- seq(0.35, 0.02, length.out=n.sets)

      text(x.coords[1], 0.37, "Underexpressed", cex=0.8, adj=0, font=3)
      text(x.coords[1], y.coords, c(1:n.genes), adj=0)
      text(x.coords[2], y.coords, round(top.gs.score, 2), cex=0.6, adj=0)
      text(x.coords[3], y.coords, format(top.gs.p, digits=1), cex=0.6, adj=0)

      text(x.coords[4], y.coords, sapply(gs.def.list[names(top.gs.score)],
                                         function(x) { length(x$Genes) }), cex=0.6, adj=0)

      text(x.coords[5], y.coords, gs.def.list.categories[names(top.gs.score)], cex=0.6, adj=0)
      text(x.coords[6], y.coords, names(top.gs.score), cex=0.6, adj=0)

      if (preferences$geneset.analysis.exact)
      {
        p <- spot.list.samples[[m]]$GSZ.p.value

        fdrtool.result <- suppressWarnings(fdrtool(p, statistic="pvalue", verbose=FALSE, plot=FALSE))
        fdr.spot.list.samples <- fdrtool.result$lfdr
        Fdr.spot.list.samples <- fdrtool.result$qval

        n.0.spot.list.samples <- fdrtool.result$param[1,"eta0"]
        perc.DE.spot.list.samples <- 1 - n.0.spot.list.samples

        par(mar=c(3,6,2,6))

        hist(p, bre=20, freq=FALSE, xlab="p-value", ylab="", main="p-values",
             las=1, cex.main=1.5, cex.lab=1, cex.axis=1)

        box()
        mtext("Density", side=2, line=3, cex=1)
        mtext("FDR", side=4, line=3, cex=1)
        mtext(paste("%DE =", round(perc.DE.spot.list.samples ,2)), line=-1.2, cex=0.5)

        abline(h=n.0.spot.list.samples , col="gray", lwd=2)

        par(new=TRUE)
        plot(0, type="n", xlim=c(0,1), ylim=c(0,1), xlab="", ylab="", axes=FALSE)
        axis(4, seq(0, 1, 0.2), seq(0, 1, 0.2), las=1, cex.axis=1)
        o = order(p)
        lines(p[o], Fdr.spot.list.samples[o], lty=2, lwd=2)
        lines(p[o], fdr.spot.list.samples[o], lty=3, lwd=3)

        legend("topright",
               c("p", expression(eta[0]), "Fdr", "fdr"),
               col=c("black","gray","black","black"),
               lty=c(1,1,2,3), lwd=c(1,1,1,2), cex=0.7)
      }
    }

    ## Spot Sheets
    if (length(spot.list.samples[[m]]$spots) <= 0)
    {
      next
    }

    for (spot.i in seq_along(spot.list.samples[[m]]$spots))
    {
      if (length(spot.list.samples[[m]]$spots[[spot.i]]$genes) <= 1)
      {
        next
      }

      spot.genes <- spot.list.samples[[m]]$spots[[spot.i]]$genes
      spot.metagenes <- spot.list.samples[[m]]$spots[[spot.i]]$metagenes

      n.genes <- min(20, length(spot.genes))
      local.p <- p.g.m[spot.genes, m]

      fdrtool.result <- suppressWarnings(fdrtool(local.p, statistic="pvalue", verbose=FALSE, plot=FALSE))

      local.fdr <- fdrtool.result$lfdr
      local.Fdr <- fdrtool.result$qval

      names(local.fdr) <- spot.genes
      names(local.Fdr) <- spot.genes

      local.n.0 <- fdrtool.result$param[1,"eta0"]
      local.perc.DE <- 1 - local.n.0

      o <- names(sort(local.p))[1:n.genes]

      layout(matrix(c(1,2,0,1,3,0,4,4,5,6,6,7), 3, 4), widths=c(1,1,2,2), heights=c(2,1,1))

      par(mar=c(0,0,0,0))
      plot(0, type="n", axes=FALSE, xlab="", ylab="", xlim=c(0,1), ylim=c(0,1))
      text(0.1, 0.94, colnames(indata)[m] , cex=3, adj=0)
      text(0.1, 0.8, "Local Summary" , cex=1.8, adj=0)
      text(0.1, 0.7,  paste("%DE =", round(local.perc.DE, 2)), adj=0)
      text(0.1, 0.65, paste("# metagenes =", length(spot.metagenes)), adj=0)
      text(0.1, 0.6, paste("# genes =", length(spot.genes)), adj=0)

      if (preferences$geneset.analysis)
      {
        text(0.1, 0.55, paste("# genes in genesets =",
                              length(intersect(unique(unlist(gs.def.list)),
                                               gene.ids[spot.genes]))), adj=0)
      }

      all.fdr.genes <- names(which(local.fdr < 0.1))
      plus.fdr.genes <- which(indata[all.fdr.genes, m] > 0)
      minus.fdr.genes <- which(indata[all.fdr.genes, m] <= 0)

      text(0.1, 0.475, paste("# genes with fdr < 0.1  =",
                             length(all.fdr.genes)," (",
                             length(plus.fdr.genes), "+ /",
                             length(minus.fdr.genes), " -)"), adj=0)

      all.fdr.genes <- names(which(local.fdr < 0.05))
      plus.fdr.genes <- which(indata[all.fdr.genes, m] > 0)
      minus.fdr.genes <- which(indata[all.fdr.genes, m] <= 0)

      text(0.1, 0.425, paste("# genes with fdr < 0.05  =",
                             length(all.fdr.genes)," (",
                             length(plus.fdr.genes), "+ /",
                             length(minus.fdr.genes), " -)"), adj=0)

      all.fdr.genes <- names(which(local.fdr < 0.01))
      plus.fdr.genes <- which(indata[all.fdr.genes, m] > 0)
      minus.fdr.genes <- which(indata[all.fdr.genes, m] <= 0)

      text(0.1, 0.375, paste("# genes with fdr < 0.01 =",
                             length(all.fdr.genes)," (",
                             length(plus.fdr.genes), "+ /",
                             length(minus.fdr.genes), " -)"), adj=0)

      text(0.1, 0.275,
           paste("<r> metagenes =", round(mean(cor(t(metadata[spot.metagenes,]))), 2)),
           adj=0)

      if (length(spot.genes) < 2000)
      {
        text(0.1, 0.225,
             paste("<r> genes =", round(mean(cor(t(indata[spot.genes,]))), 2)),
             adj=0)
      }

      text(0.1, 0.15, paste("<FC> =", round(mean(indata[spot.genes,m]), 2)), adj=0)
      text(0.1, 0.1,  paste("<t-score> =", round(mean(t.g.m[spot.genes,m]), 2)), adj=0)
      text(0.1, 0.05, paste("<p-value> =", round(10 ^ mean(log10(p.g.m[spot.genes,m])), 2)), adj=0)
      text(0.1, 0.0,  paste("<fdr> =", round(mean(fdr.g.m[spot.genes,m]), 2)), adj=0)

      par(mar=c(2,3,3,1))

      image(matrix(metadata[,m], preferences$dim.1stLvlSom, preferences$dim.1stLvlSom),
            axes=FALSE, col = colramp(1000), main="Profile", cex.main=1.5)

      axis(1, seq(0, 1, length.out = preferences$dim.1stLvlSom/10+1),
           c(1, seq(10, preferences$dim.1stLvlSom, length.out=preferences$dim.1stLvlSom/10)),
           cex.axis=1.0)

      axis(2, seq(0, 1, length.out = preferences$dim.1stLvlSom/10+1),
           c(1, seq(10, preferences$dim.1stLvlSom, length.out=preferences$dim.1stLvlSom/10)),
           cex.axis=1.0, las=1)

      box()

      image(matrix(metadata[,m], preferences$dim.1stLvlSom, preferences$dim.1stLvlSom),
            axes=FALSE, col = colramp(1000), main="Spot", cex.main=1.5)

      par(new=TRUE)
      mask <- spot.list.samples[[m]]$spots[[spot.i]]$mask
      mask[which(is.na(spot.list.samples[[m]]$spots[[spot.i]]$mask))] <- 1
      mask[which(!is.na(spot.list.samples[[m]]$spots[[spot.i]]$mask))] <- NA
      image(matrix(mask, preferences$dim.1stLvlSom, preferences$dim.1stLvlSom), axes=FALSE, col="white")

      axis(1, seq(0, 1, length.out=preferences$dim.1stLvlSom/10+1),
           c(1, seq(10, preferences$dim.1stLvlSom, length.out=preferences$dim.1stLvlSom/10)),
           cex.axis=1.0)

      axis(2, seq(0, 1, length.out=preferences$dim.1stLvlSom/10+1),
           c(1, seq(10, preferences$dim.1stLvlSom, length.out=preferences$dim.1stLvlSom/10)),
           cex.axis=1.0, las=1)

      box()

      par(mar=c(0,0,0,0))

      n.genes <- min(20, length(spot.genes))
      x.coords <- c(0, 0.06, 0.2, 0.28, 0.36, 0.44, 0.52)
      y.coords <- seq(0.75, 0.02, length.out=n.genes)
      plot(0, type="n", axes=FALSE, xlab="", ylab="", xlim=c(0,1), ylim=c(0,1))

      text(0, 0.88, "Local Genelist", cex=1.8, adj=0)

      text(x.coords, rep(c(0.82, 0.80), 4)[1:7],
           c("Rank", "ID", "log(FC)", "p-value", "fdr", "Metagene", "Description"),
           cex=1, adj=0)

      text(x.coords[1], y.coords, c(1:n.genes), adj=0)
      text(x.coords[2], y.coords, o, cex=0.6, adj=0)
      rect(x.coords[3]-0.02, y.coords[1]+0.01, 1, 0, border="white", col="white")
      text(x.coords[3], y.coords, round(indata[o, m], 2), cex=0.6, adj=0)
      text(x.coords[4], y.coords, format(local.p[o], digits=1), cex=0.6, adj=0)
      text(x.coords[5], y.coords, format(local.fdr[o], digits=1), cex=0.6, adj=0)
      text(x.coords[6], y.coords, gene.coordinates[o], cex=0.6, adj=0)
      text(x.coords[7], y.coords, gene.descriptions[o], cex=0.6, adj=0)

      par(mar=c(3,6,2,6))

      hist(local.p, bre=20, freq=FALSE, xlab="p-value", ylab="", main="p-values",
           las=1, cex.main=1.5, cex.lab=1, cex.axis=1)

      box()
      mtext("Density", side=2, line=3, cex=1)
      mtext("FDR", side=4, line=3, cex=1)
      mtext(paste("%DE =", round(local.perc.DE, 2)), line=-1.2, cex=0.5)

      abline(h = local.n.0 , col="gray", lwd=2)

      par(new=TRUE)
      plot(0, type="n", xlim=c(0,1), ylim=c(0,1), xlab="", ylab="", axes=FALSE)
      axis(4, seq(0, 1, 0.2), seq(0, 1, 0.2), las=1, cex.axis=1)

      o <- order(local.p)
      lines(local.p[o], local.Fdr[o], lty=2, lwd=2)
      lines(local.p[o], local.fdr[o], lty=3, lwd=3)

      legend("topright",
             c("p", expression(eta[0]), "Fdr", "fdr"),
             col=c("black","gray","black","black"), lty=c(1,1,2,3),
             lwd=c(1,1,1,2), cex=0.7)

      if (preferences$geneset.analysis && preferences$geneset.analysis.samplespots)
      {
        par(mar=c(0,0,0,0))

        n.sets <- 40
        x.coords <- c(0, 0.1, 0.18, 0.30, 0.39, 0.47)
        y.coords <- seq(0.75, 0.02, length.out=n.sets)

        if (spot.list.samples[[m]]$spots[[spot.i]]$type == "overexpressed")
        {
          top.gs.score <-
            sort(spot.list.samples[[m]]$spots[[spot.i]]$GSZ.score, decreasing=TRUE)[1:n.sets]

          top.gs.p <-
            spot.list.samples[[m]]$spots[[spot.i]]$GSZ.p.value[names(top.gs.score)]

          plot(0, type="n", axes=FALSE, xlab="", ylab="", xlim=c(0,1), ylim=c(0,1))

          text(0, 0.88, "Local Geneset Analysis", cex=1.8, adj=0)
          text(x.coords[1], 0.84, "Overexpression", cex=0.9, adj=0, font=3)
          text(x.coords, 0.79, c("Rank", "GSZ", "p-value", "#in/all", "Geneset", ""), cex=1, adj=0)

          text(x.coords[1], y.coords, c(1:n.sets), adj=0)
          text(x.coords[2], y.coords, round(top.gs.score, 2), cex=0.6, adj=0)
          text(x.coords[3], y.coords, format(top.gs.p, digits=1), cex=0.6, adj=0)

          text(x.coords[4],
               y.coords,
               paste(sapply(gs.def.list[names(top.gs.score)],
                            function(x) { length(intersect(x$Genes, gene.ids[spot.genes])) }),
                     "/",
                     sapply(gs.def.list[names(top.gs.score)],
                            function(x) { length(x$Genes)  })), cex=0.6, adj=0)

          text(x.coords[5], y.coords, gs.def.list.categories[names(top.gs.score)], cex=0.6, adj=0)
          text(x.coords[6], y.coords, names(top.gs.score), cex=0.6, adj=0)
        } else
        {
          top.gs.score <-
            sort(spot.list.samples[[m]]$spots[[spot.i]]$GSZ.score, decreasing=FALSE)[1:n.sets]

          top.gs.p <-
            spot.list.samples[[m]]$spots[[spot.i]]$GSZ.p.value[names(top.gs.score)]


          plot(0, type="n", axes=FALSE, xlab="", ylab="", xlim=c(0,1), ylim=c(0,1))

          text(0, 0.88, "Local Geneset Analysis", cex=1.8, adj=0)
          text(x.coords[1], 0.84, "Underexpression", cex=0.9, adj=0, font=3)
          text(x.coords, 0.79, c("Rank", "GSZ", "p-value", "#in/all", "Geneset", ""), cex=1, adj=0)

          text(x.coords[1], y.coords, c(1:n.sets), adj=0)
          text(x.coords[2], y.coords, round(top.gs.score, 2), cex=0.6, adj=0)
          text(x.coords[3], y.coords, format(top.gs.p, digits=1), cex=0.6, adj=0)

          text(x.coords[4],
               y.coords,
               paste(sapply(gs.def.list[names(top.gs.score)],
                            function(x) { length(intersect(x$Genes, gene.ids[spot.genes])) }),
                     "/",
                     sapply(gs.def.list[names(top.gs.score)],
                            function(x) { length(x$Genes) })), cex=0.6, adj=0)

          text(x.coords[5], y.coords, gs.def.list.categories[names(top.gs.score)], cex=0.6, adj=0)
          text(x.coords[6], y.coords, names(top.gs.score), cex=0.6, adj=0)
        }

        if (preferences$geneset.analysis.exact)
        {
          p <- spot.list.samples[[m]]$spots[[spot.i]]$GSZ.p.value

          fdrtool.result <- suppressWarnings(fdrtool(p, statistic="pvalue", verbose=FALSE, plot=FALSE))
          fdr.spot.list.samples <- fdrtool.result$lfdr
          Fdr.spot.list.samples <- fdrtool.result$qval

          n.0.spot.list.samples <- fdrtool.result$param[1,"eta0"]
          perc.DE.spot.list.samples <- 1 - n.0.spot.list.samples

          par(mar=c(3,6,2,6))

          hist(p, bre=20, freq=FALSE, xlab="p-value", ylab="", main="p-values",
               las=1, cex.main=1.5, cex.lab=1, cex.axis=1)

          box()
          mtext("Density", side=2, line=3, cex=1)
          mtext("FDR", side=4, line=3, cex=1)
          mtext(paste("%DE =", round(perc.DE.spot.list.samples ,2)), line=-1.2, cex=0.5)

          abline(h = n.0.spot.list.samples , col="gray", lwd=2)

          par(new=TRUE)
          plot(0, type="n", xlim=c(0,1), ylim=c(0,1), xlab="", ylab="", axes=FALSE)
          axis(4, seq(0, 1, 0.2), seq(0, 1, 0.2), las=1, cex.axis=1)
          o = order(p)
          lines(p[o], Fdr.spot.list.samples[o], lty=2, lwd=2)
          lines(p[o], fdr.spot.list.samples[o], lty=3, lwd=3)

          legend("topright", c("p", expression(eta[0]), "Fdr", "fdr"),
                 col=c("black","gray","black","black"), lty=c(1,1,2,3),
                 lwd=c(1,1,1,2), cex=0.7)
        }
      }
    }

    dev.off()
  }
}
