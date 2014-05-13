pipeline.groupSpecificGenesets <- function()
{
  gr <- 1

  for (gr in 1:length(unique(group.labels)))
  {

    samples.GSZ.scores <- do.call(cbind, lapply(spot.list.samples, function(x)
    {
      x$GSZ.score[names(gs.def.list)]
    }))

    gs.p.values <- apply(samples.GSZ.scores, 1, function(x)
    {
      wilcox.test(x[which(group.labels!=unique(group.labels)[gr])],
                  x[which(group.labels==unique(group.labels)[gr])],
                  alternative="less")$p.value
    })

    gs.p.values <- sort(gs.p.values)
    top.gs <- gs.p.values[1:20]

    pdf(paste(files.name,
              " - Results/Summary Sheets - Groups/Geneset Analysis/Specific GS ",
              make.names(unique(group.labels)[gr]),".pdf", sep=""), 21/2.54, 29.7/2.54)

    layout(matrix(c(1:8),4, byrow=T), widths=c(3,1))

    for (i in 1:20)
    {
      ylim <- c(-15, 20)
      par(mar=c(3,5,2,1))

      barplot(samples.GSZ.scores[names(top.gs)[i],], beside=T, las=2,
              cex.names=1.2, col=group.colors, cex.main=1, ylim=ylim,
              border=if (ncol(indata) < 80) "black" else NA,
              names.arg=if (ncol(indata)<100) colnames(indata) else rep("",ncol(indata)),
              cex.axis=1.4)

      abline(h=0,lty=2)
      title(main=names(top.gs)[i], cex.main=1.2, line=1)
      mtext("GSZ", side=2, cex=1.2, line=3)


      n.map <- matrix(0,preferences$dim.1stLvlSom,preferences$dim.1stLvlSom)
      gs.nodes <- som.nodes[names(gene.ids)[which(gene.ids %in% gs.def.list[[names(top.gs)[i]]]$Genes)]]
      n.map[as.numeric(names(table(gs.nodes)))] <- table(gs.nodes)
      n.map[which(n.map==0)] <- NA
      n.map <- matrix(n.map, preferences$dim.1stLvlSom)

      par(mar=c(5,1,3,1))

      lim <- c(1,preferences$dim.1stLvlSom) + preferences$dim.1stLvlSom*0.01*c(-1,1)

      plot(which(!is.na(n.map), arr.ind=T), xlim=lim, ylim=lim, pch=16,
            axes=F, xlab="",ylab="", xaxs="i", yaxs="i",
            cex=0.5 + na.omit(as.vector(n.map)) / max(n.map,na.rm=T) * 2.8,
            col=colramp(1000)[(na.omit(as.vector(n.map)) - min(n.map,na.rm=T)) /
                              max(1, (max(n.map,na.rm=T) - min(n.map,na.rm=T))) *
                              999 + 1])

      title(sub=paste("# features =", length(gs.def.list[[names(top.gs)[i]]]$Genes),
                      ", max =", max(n.map,na.rm=T)),line=0.5, cex.sub=1.4)

      box()
    }

    dev.off()

    fdr.res <- fdrtool(gs.p.values,statistic="pvalue",plot=F,verbose=F)
    out <- cbind(names(gs.p.values), gs.p.values, fdr.res$lfdr)
    colnames(out) = c("gene set","p-value","fdr")

    write.csv2(out, paste(files.name,
                          " - Results/Summary Sheets - Groups/Geneset Analysis/Specific GS ",
                          make.names(unique(group.labels)[gr]),".csv", sep=""), row.names=F)
  }
}
