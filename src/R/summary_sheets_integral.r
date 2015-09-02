pipeline.summarySheetsIntegral <- function()
{
  #### Summary Sheets ####
  plot.set.list <- function(set.list, main, path)
  {
    pdf(path, 29.7/2.54, 21/2.54)

    if (main != "Correlation Cluster" && main != "K-Means Cluster")
    {
      layout(matrix(c(1, 2), 1, 2), c(2, 1), 1)
      par(mar=c(5, 4, 4, 1))

      col <- if(main!="D-Clusters") colramp(1000) else colorRampPalette(c("blue2","white","red2"))(1000)
      image(matrix(set.list$overview.map, preferences$dim.1stLvlSom, preferences$dim.1stLvlSom),
            axes=FALSE, col=col, main=main, cex.main=1.5)

      mtext("landscape", 3)
      box()
      par(mar=c(5, 1, 4, 2))

      plot(0, type="n", axes=FALSE, xlab="", ylab="", xlim=c(0,1), ylim=c(0,1),
           xaxs="i", yaxs="i")

      box()
    }

    # Spot Altas "beta-colored"
    layout(matrix(c(1, 2), 1, 2), c(2, 1), 1)
    par(mar=c(5, 4, 4, 1))

    beta.scores <- sapply(set.list$spots, function(x) { x$beta.statistic$beta.score })

    image(matrix(set.list$overview.mask, preferences$dim.1stLvlSom), axes=FALSE,
          col=colramp(1000)[1 + 999 * beta.scores / max(beta.scores)],
          main=main, cex.main=1.5)

    mtext("beta-scores", 3)
    box()
    par(mar=c(5, 1, 4, 2))

    plot(0, type="n", axes=FALSE, xlab="", ylab="", xlim=c(0,1), ylim=c(0,1),
         xaxs="i", yaxs="i")

    box()
    par(new=TRUE)
    par(mar=c(12, 9, 12, 9))
    image(matrix(1:100, 1, 100), col = colramp(1000), axes=FALSE)
    axis(2, c(0,1), c(0, round(max(beta.scores),1)), las=2)
    box()

    # Spot Altas "politically"
    layout(matrix(c(1, 2), 1, 2), c(2, 1), 1)
    par(mar=c(5, 4, 4, 1))

    image(x=c(1:preferences$dim.1stLvlSom),
          y=c(1:preferences$dim.1stLvlSom),
          z=matrix(set.list$overview.mask, preferences$dim.1stLvlSom, preferences$dim.1stLvlSom),
          col=colramp(max(set.list$overview.mask, na.rm=TRUE)),
          axes=TRUE, main=main, cex.main=1.5, xlab="", ylab="", las=1)

    mtext("annotation", 3)
    box()
    par(new=TRUE)

    plot(0, type="n", axes=FALSE, xlab="", ylab="", xlim=c(0,preferences$dim.1stLvlSom),
         ylim=c(0,preferences$dim.1stLvlSom), xaxs="i", yaxs="i")

    points(do.call(rbind, lapply(set.list$spots, function(x) { x$position })),
           pch=16, cex=3, col="black")

    points(do.call(rbind, lapply(set.list$spots, function(x) { x$position })),
           pch=1, cex=3, col="white")

    text(do.call(rbind, lapply(set.list$spots, function(x) { x$position })),
         names(set.list$spots), col="white")

    par(mar=c(5, 1, 4, 2))

    plot(0, type="n", axes=FALSE, xlab="", ylab="", xlim=c(0,1), ylim=c(0,1),
         xaxs="i", yaxs="i")

    box()

    if (preferences$geneset.analysis)
    {
      top.GS <- lapply(set.list$spots, function(x) names(head(x$Fisher.p , 3)))

      leg.col <- colramp(length(set.list$spots))
      leg.col <- as.vector(sapply(leg.col, c, NA, NA))
      leg.num <- names(set.list$spots)
      leg.num <- as.vector(sapply(leg.num, c, NA, NA))

      legend(x=0.05, y=1, unlist(top.GS), cex=0.7, col = leg.col, pch=15,
             pt.cex=1.5, bty="n")

      legend(x=-0.04, y=1, legend=leg.num, cex=0.7, bty="n")
    }  else
    {
      top.GS <- rep("", 3 * length(set.list$spots))
    }

    # Spot - Sample - Heatmap
    sample.spot.expression <- matrix(NA, 0, ncol(indata))

    for (m in seq_along(set.list$spots))
    {
      mean.FC <- apply(metadata, 2, function(x)
      {
        max(x[set.list$spots[[m]]$metagenes])
      })

      sample.spot.expression <- rbind(sample.spot.expression, mean.FC)
    }

    rownames(sample.spot.expression) <- names(set.list$spots)

    sample.spot.expression.image <-
      if (nrow(sample.spot.expression) > 1)
      {
        t(sample.spot.expression[nrow(sample.spot.expression):1,])
      } else
      {
        as.matrix(sample.spot.expression[nrow(sample.spot.expression):1,])
      }

    layout(matrix(c(0,2,0,3,1,0,0,4,5), 3, 3), heights=c(0.8,6,2), widths=c(0.5,5,3))

    par(mar=c(0,0,0,0))

    image(1:ncol(indata),
          1:nrow(sample.spot.expression),
          sample.spot.expression.image,
          col=colorRampPalette(c("blue4","blue","gray90","orange","red4"))(1000),
          axes=FALSE, ylim=0.5+c(0,nrow(sample.spot.expression)), yaxs="i", xlab="", ylab="",
          zlim=max(max(sample.spot.expression.image),-min(sample.spot.expression.image))*c(-1,1))

    box()

    if (ncol(indata)<100)
    {
      axis(1, 1:ncol(indata), labels=colnames(indata), las=2, line=-0.5, tick=0, cex.axis=1.4)
    }

    plot(0, type="n", xlab="", ylab="", axes=FALSE, xlim=c(0,1),
         ylim=0.5+c(0,nrow(sample.spot.expression)), yaxs="i")

    text(0.7, nrow(sample.spot.expression):1, rownames(sample.spot.expression),
         adj=1, cex=1.8)

    par(mar=c(1,0,2,0))

    if (length(unique(group.labels)) > 1)
    {
      image(cbind(1:ncol(indata)), col = group.colors, axes = FALSE)
      box()
    } else
    {
      frame()
    }

    par(mar=c(0,0,0,0))

    plot(0, type="n", xlab="", ylab="", axes=FALSE, xlim=c(0,1),
         ylim=0.5+c(0,nrow(sample.spot.expression)), yaxs="i")

    pos <- as.vector(sapply(c(1:nrow(sample.spot.expression)),
                            function(x) { c(x-0.26, x, x+0.26) }))

    text(0.05,
         rev(pos),
         unlist(lapply(set.list$spots, function(x) { names(head(x$Fisher.p , 3)) })),
         adj=0, cex=1) #0.6

    par(mar=c(5,2,4,2))

    image(matrix(1:100, 100, 1),
          col=colorRampPalette(c("blue4","blue","gray90","orange","red4"))(1000),
          axes=FALSE, xlab="")

    axis(1, round(max(max(sample.spot.expression.image),
                      -min(sample.spot.expression.image)) * c(-1,1), 1),
         at=c(0,1), las=2, tick=FALSE, pos=-0.8, cex.axis=1.4)

    mtext(expression(paste("<",Delta,"e", '' ^ meta, ">")), side=1, line=0.5)
    box()

    # Spot - SampleGSZ - Heatmap
    if (preferences$geneset.analysis)
    {
      sample.spot.GSZ <- sapply(spot.list.samples, function(x)
      {
        return(x$GSZ.score[unlist(top.GS)])
      })

      sample.spot.GSZ.image <-
        if (nrow(sample.spot.GSZ) > 1)
        {
          t(sample.spot.GSZ[nrow(sample.spot.GSZ):1,])
        } else
        {
          as.matrix(sample.spot.GSZ[nrow(sample.spot.GSZ):1,])
        }

      sample.spot.GSZ.image[which(is.na(sample.spot.GSZ.image))] <- 0
      layout(matrix(c(0,2,0,3,1,0,0,4,5), 3, 3), heights=c(0.8,6,2), widths=c(0.5,5,3))
      par(mar=c(0,0,0,0))

      image(1:ncol(indata), 1:nrow(sample.spot.GSZ), sample.spot.GSZ.image,
            axes=FALSE, ylim=0.5+c(0,nrow(sample.spot.GSZ)), yaxs="i", xlab="", ylab="",
            col=colorRampPalette(c("blue4","blue","gray90","orange","red4"))(1000),
            zlim=max(max(sample.spot.GSZ.image),-min(sample.spot.GSZ.image))*c(-1,1))

      box()

      if (ncol(indata)<100)
      {
        axis(1, 1:ncol(indata), labels=colnames(indata), las=2, line=-0.5,
             tick=0, cex.axis=1.4)
      }

      plot(0, type="n", xlab="", ylab="", axes=FALSE, xlim=c(0,1),
           ylim=0.5+c(0,nrow(sample.spot.expression)), yaxs="i")

      text(0.7, nrow(sample.spot.expression):1, names(set.list$spots),
           adj=1, cex=1.8)

      par(mar=c(1,0,2,0))

      if (length(unique(group.labels)) > 1)
      {
        image(cbind(1:ncol(indata)), col=group.colors, axes=FALSE)
        box()
      } else
      {
        plot(0, type="n", xlab="", ylab="", axes=FALSE)
      }

      par(mar=c(0,0,0,0))

      plot(0, type="n", xlab="", ylab="", axes=FALSE, xlim=c(0,1),
           ylim=0.5+c(0,nrow(sample.spot.GSZ)), yaxs="i")

      pos <- as.vector(sapply(c(1:nrow(sample.spot.GSZ)),
                              function(x) { c(x-0.26, x, x+0.26) }))

      text(0.05, rev(c(1:nrow(sample.spot.GSZ)) + c(0.16,0,-0.16)),
           unlist(lapply(set.list$spots, function(x) { names(head(x$Fisher.p , 3)) })),
           adj=0, cex=1) #0.6

      par(mar=c(5,2,4,2))

      image(matrix(1:100, 100, 1),
            col=colorRampPalette(c("blue4","blue","gray90","orange","red4"))(1000),
            axes=FALSE, xlab="")

      axis(1,
           round(max(max(sample.spot.GSZ.image),-min(sample.spot.GSZ.image))*c(-1,1), 1),
           at=c(0,1), las=2, tick=FALSE, pos=-0.8, cex.axis=1.4)

      mtext("sample GSZ", side=1, line=0.5)
      box()
    }

    # Individual spot sheets
    for (m in seq_along(set.list$spots))
    {
      if (main %in% c("Sample-Underexpression"))
      {
        sample.with.spot <- set.list$spotdata[m,] < -sd(as.vector(set.list$spotdata))
      }  else
      {
        sample.with.spot <- set.list$spotdata[m,] > sd(as.vector(set.list$spotdata))
      }

      layout(matrix(c(1,2,4,1,3,4,5,5,6,7,7,8), 3, 4), widths=c(1,1,2,2), heights=c(2,1,1))

      par(mar=c(0,0,0,0))
      plot(0, type="n", axes=FALSE, xlab="", ylab="", xlim=c(0,1), ylim=c(0,1))
      text(0.1, 0.94, main , cex=2.6, adj=0)
      text(0.1, 0.8, paste("Spot Summary:", names(set.list$spots)[m]) , cex=1.8, adj=0)
      text(0.1, 0.7, paste("# metagenes =", length(set.list$spots[[m]]$metagenes)), adj=0)
      text(0.1, 0.66, paste("# genes =", length(set.list$spots[[m]]$genes)), adj=0)

      text(0.1, 0.55,
           paste("<r> metagenes =",
                 round(mean(cor(t(metadata[set.list$spots[[m]]$metagenes,]))), 2)), adj=0)

      if (length(set.list$spots[[m]]$genes) < 1000)
      {
        suppressWarnings({
          try({
            text(0.1, 0.51,
                 paste("<r> genes =",
                       round(mean(cor(t(indata[set.list$spots[[m]]$genes,]))), 2)),
                 adj=0)
          }, silent=TRUE)
        })
      }

      text(0.1, 0.47,
           paste("beta: r2=",
                 round(set.list$spots[[m]]$beta.statistic$beta.score,2),
                 " /  log p=",
                 round(log10(set.list$spots[[m]]$beta.statistic$beta.significance), 2)),
           adj=0)

      text(0.1, 0.39,
           paste("# samples with spot =",
                 sum(sample.with.spot), "(",
                 round(100 * sum(sample.with.spot) / ncol(indata), 1), "%)"), adj=0)

      if (length(unique(group.labels)) > 1 && sum(sample.with.spot) > 0)
      {
        group.table <- table(group.labels[sample.with.spot])[unique(group.labels)]
        group.table <- group.table[which(!is.na(group.table))]

        for (g in seq_along(group.table))
        {
          text(0.15, 0.39-g*0.04,
               paste(names(group.table)[g], ":", group.table[g], "(",
                     round(100 * group.table[g]/sum(group.labels == names(group.table)[g]), 1),
                     "%)"), adj=0,  col=group.colors[match(names(group.table)[g], group.labels)])
         }
      }

      par(mar=c(2,3,3,1))

      col <- if(main!="D-Clusters") colramp(1000) else colorRampPalette(c("blue2","white","red2"))(1000)
      image(matrix(set.list$overview.map, preferences$dim.1stLvlSom, preferences$dim.1stLvlSom),
            axes=FALSE, col=col, main="Overview Map", cex.main=1.5)

      axis(1, seq(0, 1, length.out = preferences$dim.1stLvlSom/10+1),
           c(1, seq(10, preferences$dim.1stLvlSom, length.out = preferences$dim.1stLvlSom/10)),
           cex.axis=1.0)

      axis(2, seq(0, 1, length.out = preferences$dim.1stLvlSom/10+1),
           c(1, seq(10, preferences$dim.1stLvlSom, length.out = preferences$dim.1stLvlSom/10)),
           cex.axis=1.0, las=1)

      box()

      image(matrix(set.list$overview.map, preferences$dim.1stLvlSom, preferences$dim.1stLvlSom),
            axes=FALSE, col=col, main="Spot", cex.main=1.5)

      par(new=TRUE)

      mask <- set.list$spots[[m]]$mask
      mask[which(is.na(set.list$spots[[m]]$mask))] <- 1
      mask[which(!is.na(set.list$spots[[m]]$mask))] <- NA

      image(matrix(mask, preferences$dim.1stLvlSom, preferences$dim.1stLvlSom),
            axes=FALSE, col = "white")

      axis(1, seq(0, 1, length.out = preferences$dim.1stLvlSom/10+1),
           c(1, seq(10, preferences$dim.1stLvlSom, length.out = preferences$dim.1stLvlSom/10)),
           cex.axis=1.0)

      axis(2, seq(0, 1, length.out = preferences$dim.1stLvlSom/10+1),
           c(1, seq(10, preferences$dim.1stLvlSom, length.out = preferences$dim.1stLvlSom/10)),
           cex.axis=1.0, las=1)

      box()

      # Spot Profile Plot
      par(mar=c(8,3,1,1))

      barplot(set.list$spotdata[m,], col=group.colors, main="",
              names.arg=if (ncol(indata)<100) colnames(indata) else rep("",ncol(indata)),
              las=2, cex.main=1, cex.lab=1, cex.axis=1, cex.names=0.8,
              border=if (ncol(indata) < 80) "black" else NA)

      box()

      if (length(set.list$spots[[m]]$genes) > 0 && length(set.list$spots[[m]]$genes) < 5000)
      {
        # Spot Genelist
        r.genes <- sapply(set.list$spots[[m]]$genes, function(x)
        {
          gene <- indata[x,]
          return(suppressWarnings(cor(gene, set.list$spotdata[m,])))
        })

        e.max <- apply(indata[set.list$spots[[m]]$genes, ,drop=FALSE], 1, max)
        e.min <- apply(indata[set.list$spots[[m]]$genes, ,drop=FALSE], 1, min)

        if (main %in% c("Sample-Underexpression","Metagene Minima"))
        {
          o <- names(sort(e.min, decreasing=FALSE))
        }  else
        {
          o <- names(sort(e.max, decreasing=TRUE))
        }

        n.genes <- 20
        o <- o[1:min(n.genes,length(o))]

        par(mar=c(0,0,0,0))

        x.coords <- c(0, 0.06, 0.2, 0.28, 0.36, 0.44, 0.52)
        y.coords <- seq(0.75, 0.02, length.out=length(o))

        plot(0, type="n", axes=FALSE, xlab="", ylab="", xlim=c(0,1), ylim=c(0,1))

        text(0, 0.88, "Spot Genelist", cex=1.8, adj=0)

        text(x.coords, rep(c(0.82, 0.80), 4)[1:7],
             c("Rank", "ID", "max e", "min e", "r", "Symbol", "Description"),
             cex=1, adj=0)

        text(x.coords[1], y.coords, c(seq_along(o)), adj=0)
        text(x.coords[2], y.coords, o, cex=0.6, adj=0)
        rect(x.coords[3]-0.02, y.coords[1]+0.01, 1, 0, border="white", col="white")
        text(x.coords[3], y.coords, round(e.max[o], 2), cex=0.6, adj=0)
        text(x.coords[4], y.coords, round(e.min[o], 2), cex=0.6, adj=0)
        text(x.coords[5], y.coords, round(r.genes[o], 2), cex=0.6, adj=0)
        text(x.coords[6], y.coords, gene.names[o], cex=0.6, adj=0)
        text(x.coords[7], y.coords, gene.descriptions[o], cex=0.6, adj=0)
      } else
      {
        frame()
      }

      plot(0, type="n", axes=FALSE, xlab="", ylab="", xlim=c(0,1), ylim=c(0,1))

      if (preferences$geneset.analysis)
      {
        n.sets <- 40
        top.gs.p <- sort(set.list$spots[[m]]$Fisher.p)[1:n.sets]
        par(mar=c(0,0,0,0))

        x.coords <- c(0, 0.1, 0.23, 0.34, 0.4)
        y.coords <- seq(0.75, 0.02, length.out=n.sets)

        plot(0, type="n", axes=FALSE, xlab="", ylab="", xlim=c(0,1), ylim=c(0,1))

        text(0, 0.88, "Geneset Overrepresentation", cex=1.8, adj=0)

        text(x.coords, 0.82, c("Rank", "p-value", "#in/all", "Geneset", ""),
             cex=1, adj=0)

        text(x.coords[1], y.coords, c(1:n.sets), adj=0)
        text(x.coords[2], y.coords, format(top.gs.p, digits=1), cex=0.6, adj=0)

        text(x.coords[3], y.coords, paste (sapply(gs.def.list[names(top.gs.p)], function(x)
        {
          length(intersect(x$Genes, gene.ids[set.list$spots[[m]]$genes]))
        }), "/", sapply(gs.def.list[names(top.gs.p)], function(x)
        {
          length(x$Genes)
        })), cex=0.6, adj=0)

        text(x.coords[4], y.coords, gs.def.list.categories[names(top.gs.p)], cex=0.6, adj=0)
        rect(x.coords[5]-0.01, y.coords+0.01, 1, 0, border="white", col="white")
        text(x.coords[5], y.coords, names(top.gs.p), cex=0.6, adj=0)

        try.res <- try({
          fdrtool.result <- suppressWarnings({
            fdrtool(set.list$spots[[m]]$Fisher.p, statistic="pvalue", verbose=FALSE, plot=FALSE)
          })
        }, silent=TRUE)

        if (class(try.res) != "try-error")
        {
          fdr.spot.list.samples <- fdrtool.result$lfdr
          Fdr.spot.list.samples <- fdrtool.result$qval

          n.0.spot.list.samples <- fdrtool.result$param[1,"eta0"]
          perc.DE.spot.list.samples <-1 - n.0.spot.list.samples

          par(mar=c(3,6,2,6))

          hist(set.list$spots[[m]]$Fisher.p, bre=20, freq=FALSE, xlab="p-value",
               ylab="", main="p-values", las=1, cex.main=1.5, cex.lab=1, cex.axis=1)

          box()
          mtext("Density", side=2, line=3, cex=1)
          mtext("FDR", side=4, line=3, cex=1)
          mtext(paste("%DE =", round(perc.DE.spot.list.samples ,2)), line=-1.2, cex=0.5)

          abline(h = n.0.spot.list.samples , col="gray", lwd=2)

          par(new=TRUE)
          plot(0, type="n", xlim=c(0,1), ylim=c(0,1), xlab="", ylab="", axes=FALSE)
          axis(4, seq(0, 1, 0.2), seq(0, 1, 0.2), las=1, cex.axis=1)
          o = order(set.list$spots[[m]]$Fisher.p)
          lines(set.list$spots[[m]]$Fisher.p[o], Fdr.spot.list.samples[o], lty=2, lwd=2)
          lines(set.list$spots[[m]]$Fisher.p[o], fdr.spot.list.samples[o], lty=3, lwd=3)

          legend("topright", c("p", expression(eta[0]), "Fdr", "fdr"),
                 col=c("black","gray","black","black"), lty=c(1,1,2,3), lwd=c(1,1,1,2), cex=0.7)
        }
        
        ## Splitted Genesets Sheet
        n.sets <- 20
        n.cat <- length(table(gs.def.list.categories))
        par(mfrow=c(ceiling(n.cat/3), min(n.cat, 3)))
        
        for (i in names(table(gs.def.list.categories)))
        {
          top.gs.p <-
            sort(set.list$spots[[m]]$Fisher.p[names(which(gs.def.list.categories == i))])[1:n.sets]
          
          x.coords <- c(0.05, 0.15, 0.28, 0.39, 0.45)
          y.coords <- seq(0.88, 0.05, length.out=n.sets)
          
          par(mar=c(0,0,0,0))
          
          plot(0, type="n", axes=FALSE, xlab="", ylab="", xlim=c(0,1),
               ylim=c(0,1), xaxs="i", yaxs="i")
          
          text(x.coords[1], 0.97, i, cex=2, adj=0)
          text(x.coords, 0.92, c("Rank", "p-value", "#in/all", "Geneset", ""), cex=1, adj=0)
          text(x.coords[1], y.coords, c(1:n.sets), adj=0)
          text(x.coords[2], y.coords, format(top.gs.p, digits=1), cex=0.6, adj=0)
          
          text(x.coords[3], y.coords,
               paste(sapply(gs.def.list[names(top.gs.p)], function(x)
               {
                 length(intersect(x$Genes, gene.ids[set.list$spots[[m]]$genes]))
               }), "/",
               sapply(gs.def.list[names(top.gs.p)], function(x)
               {
                 length(x$Genes)
               })), cex=0.6, adj=0)
          
          text(x.coords[4], y.coords, names(top.gs.p), cex=0.6, adj=0)
        }
        
      }
    }

    dev.off()
  }

  csv.set.list <- function(set.list, main, path)
  {
    for (m in seq_along(set.list$spots))
    {
      basename <- paste(main, " ", names(set.list$spots)[m], ".csv", sep="")

      if (length(set.list$spots[[m]]$genes) <= 0 || length(set.list$spots[[m]]$genes) > 2000)
      {
        next
      }

      ## CSV Table
      r.genes <- sapply(set.list$spots[[m]]$genes, function(x)
      {
        gene <- indata[x,]
        return(suppressWarnings(cor(gene, set.list$spotdata[m,])))
      })

      r.t <- r.genes / sqrt((1-r.genes^2) / (ncol(indata)-2))
      r.p <- 1 - pt(r.t, ncol(indata)-2)


      e.max <- apply(indata[set.list$spots[[m]]$genes, ,drop=FALSE], 1, max)
      e.min <- apply(indata[set.list$spots[[m]]$genes, ,drop=FALSE], 1, min)

      if (main %in% c("Sample-Underexpression","Metagene Minima"))
      {
        o <- names(sort(e.min, decreasing=FALSE))
      }  else
      {
        o <- names(sort(e.max, decreasing=TRUE))
      }

      out <- data.frame(Rank=c(seq_along(set.list$spots[[m]]$genes)),
                        ID=o,
                        Symbol=gene.names[o])

      if (any(grep("GenesetSOM",files.name)))
      {
        out <-cbind(out, Type = gs.def.list.categories[o])
      }

      high.low.threshold <- mean(indata.gene.mean)

      out <- cbind(out,
                   "mean expression"=indata.gene.mean[o],
                   "SD"=apply(indata[o, ,drop=FALSE], 1, sd),
                   "max delta e"=apply(indata[o, ,drop=FALSE], 1, max),
                   "min delta e"=apply(indata[o, ,drop=FALSE], 1, min),
                   "% high expressed"=round(apply(indata[o, ,drop=FALSE] + indata.gene.mean[o], 1, function(x)  sum(x > high.low.threshold) / length(x) * 100)),
                   "correlation"=r.genes[o],
                   "->t.score"=r.t[o],
                   "->p.value"=r.p[o],
                   "Metagene"=gene.coordinates[o],
                   "Chromosome"=gene.positions[o],
                   "Description"=gene.descriptions[o])

      write.csv2(out, file.path(path, basename))
    }
  }

  # directories to store the results
  dirnames <- c("pdf"=as.character(output.paths["Summary Sheets Integral"]),
                "csv"=file.path(output.paths["CSV"], "Spot Lists"))

  for (dirname in dirnames)
  {
    dir.create(dirname, showWarnings=FALSE, recursive=TRUE)
  }
  
  
  # pdf sheets to generate
  pdf.sheets <- list(
    list("Overexpression.pdf", "Overexpression Spots", spot.list.overexpression),
    list("Underexpression.pdf", "Underexpression Spots", spot.list.underexpression),
    list("Correlation Cluster.pdf", "Correlation Clusters", spot.list.correlation),
    list("K-Means Cluster.pdf", "K-Means Clusters", spot.list.kmeans),
    list("D-Clusters.pdf", "D-Clusters", spot.list.dmap))

  # csv sheets to generate
  csv.sheets <- list(
    list("Overexpression Spots", spot.list.overexpression),
    list("Underexpression Spots", spot.list.underexpression),
    list("Correlation Clusters", spot.list.correlation),
    list("K-Means Clusters", spot.list.kmeans),
    list("D-Clusters", spot.list.dmap))

  # generate group expression sheets?
  if (length(unique(group.labels)) > 1)
  {
    pdf.sheets[[length(pdf.sheets)+1]] <-
      list("Group Overexpression.pdf", "Group Overexpression Spots", spot.list.group.overexpression)

    csv.sheets[[length(csv.sheets)+1]] <-
      list("Group Overexpression Spots", spot.list.group.overexpression)
  }

  pdf.sheets <- lapply(pdf.sheets, function(x)
  {
    list(fn=plot.set.list,
         args=list(path=file.path(dirnames[["pdf"]], x[[1]]), main=x[[2]], set.list=x[[3]]))
  })

  csv.sheets <- lapply(csv.sheets, function(x)
  {
    list(fn=csv.set.list,
         args=list(path=dirnames[["csv"]], main=x[[1]], set.list=x[[2]]))
  })




  util.info("Writing:", file.path(dirnames[["pdf"]], "*.pdf"))

  dummy = sapply( seq_along(pdf.sheets), function(i, pdf.sheets)
  {
    do.call(pdf.sheets[[i]]$fn, pdf.sheets[[i]]$args)
  }, pdf.sheets)




  util.info("Writing:", file.path(dirnames[["csv"]], "*.csv"))
  
  dummy = sapply( seq_along(csv.sheets), function(i, csv.sheets)
  {
    do.call(csv.sheets[[i]]$fn, csv.sheets[[i]]$args)
  }, csv.sheets)



}
