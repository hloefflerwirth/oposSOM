pipeline.geneLists <- function()
{
  dirnames <- c("global"=file.path(output.paths["CSV"], "Gene Lists - Global"),
                "local"=file.path(output.paths["CSV"], "Gene Lists - Local"),
                "set"=file.path(output.paths["CSV"], "Gene Set Lists - Global"))

  for (dirname in dirnames)
  {
    dir.create(dirname, showWarnings=F)
  }

  #### Global Gene Lists ####
  util.info("Writing:", file.path(dirnames["global"], "*.csv"))

  for (m in 1:ncol(indata))
  {
    o <- order(p.g.m[,m])

    out <- data.frame(Rank=c(1:nrow(indata)),
                      ID=rownames(indata)[o],
                      Symbol=gene.names[o])

    if (any(grep("GenesetSOM", files.name)))
    {
      out <- cbind(out, Type=sapply(gs.def.list[rownames(indata)[o]],
                                    function(x) { x$Type }))
    }

    out <- cbind(out,
                 logFC=indata[o, m],
                 WAD=WAD.g.m[o, m],
                 T.Score=t.g.m[o, m],
                 p.value=p.g.m[o, m],
                 fdr=fdr.g.m[o, m],
                 Fdr=Fdr.g.m[o, m],
                 Metagene=genes.coordinates[o],
                 Chromosome=gene.positions[rownames(indata)[o]],
                 Description=gene.descriptions[o])

    basename <- paste(colnames(indata)[m], ".csv", sep="")
    f <- file(file.path(dirnames["global"], basename), "w")

    writeLines("Sample Summary:", f)
    writeLines("", f)
    writeLines(paste("%DE:" ,"", round(perc.DE.m[colnames(indata)[m]], 2), sep=";"), f)
    writeLines(paste("#genes with fdr < 0.2" ,"", length(which(fdr.g.m[,m] < 0.2)) , sep=";"), f)
    writeLines(paste("#genes with fdr < 0.1" ,"", length(which(fdr.g.m[,m] < 0.1)) , sep=";"), f)
    writeLines(paste("#genes with fdr < 0.05" ,"", length(which(fdr.g.m[,m] < 0.05)) , sep=";"), f)
    writeLines(paste("#genes with fdr < 0.01" ,"", length(which(fdr.g.m[,m] < 0.01)) , sep=";"), f)
    writeLines("", f)
    writeLines(paste("<FC> =", round(mean(indata[,m]), 2))  ,f)
    writeLines(paste("<t-score> =", round(mean(t.g.m[,m]), 2))  ,f)
    writeLines(paste("<p-value> =", round(10 ^ mean(log10(p.g.m[,m])), 2))  ,f)
    writeLines(paste("<fdr> =", round(mean(fdr.g.m[,m]), 2))  ,f)
    writeLines("", f);  writeLines("", f);  writeLines("", f)
    writeLines("Gene Statistics", f)
    writeLines("", f)
    write.csv2(out, file=f, row.names=F)

    close(f)
  }

  #### Local Gene Lists ####
  util.info("Writing:", file.path(dirnames["local"], "*.csv"))

  for (m in 1:ncol(indata))
  {
    if (length(GS.infos.samples[[m]]$spots) <= 0)
    {
      next
    }

    for (spot.i in 1:length(GS.infos.samples[[m]]$spots))
    {
      if (length(GS.infos.samples[[m]]$spots[[spot.i]]$genes) <= 1)
      {
        next
      }

      spot.genes <- GS.infos.samples[[m]]$spots[[spot.i]]$genes
      spot.metagenes <- GS.infos.samples[[m]]$spots[[spot.i]]$metagenes

      p <- p.g.m[spot.genes, m]
      fdrtool.result <- suppressWarnings(fdrtool(p, statistic="pvalue", verbose=F, plot=F))
      fdr.spot <- fdrtool.result$lfdr
      Fdr.spot <- fdrtool.result$qval

      names(fdr.spot) <- spot.genes
      names(Fdr.spot) <- spot.genes

      n.0.spot <- fdrtool.result$param[1,"eta0"]
      perc.DE.spot <- 1 - n.0.spot

      o <- names(sort(p))

      out <- data.frame(Rank=c(1:length(spot.genes)),
                        ID=o,
                        Symbol=gene.names[o])

      if (any(grep("GenesetSOM", files.name)))
      {
        out <- cbind(out, Type=sapply(gs.def.list[o], function(x) { x$Type }))
      }

      out <- cbind(out,
                   logFC=indata[o, m],
                   WAD=WAD.g.m[o, m],
                   T.Score=t.g.m[o, m],
                   p.value=p.g.m[o, m],
                   fdr=fdr.spot[o],
                   Fdr=Fdr.spot[o],
                   Metagene=genes.coordinates[o],
                   Chromosome=gene.positions[o],
                   Description=gene.descriptions[o])

      basename <- paste(colnames(indata)[m], ".", spot.i, ".csv", sep="")
      f <- file(file.path(dirnames["local"], basename), "w")

      writeLines("Spot Summary:", f)
      writeLines("", f)
      writeLines(paste("%DE:" ,"", round(perc.DE.spot,2), sep=";"), f)
      writeLines(paste("# metagenes:" ,"", length(spot.metagenes), sep=";"), f)
      writeLines(paste("# genes:" ,"", length(spot.genes), sep=";"), f)
      writeLines(paste("#genes with fdr < 0.2" ,"", length(which(fdr.spot < 0.2)) , sep=";"), f)
      writeLines(paste("#genes with fdr < 0.1" ,"", length(which(fdr.spot < 0.1)) , sep=";"), f)
      writeLines(paste("#genes with fdr < 0.05" ,"", length(which(fdr.spot < 0.05)) , sep=";"), f)
      writeLines(paste("#genes with fdr < 0.01" ,"", length(which(fdr.spot < 0.01)) , sep=";"), f)
      writeLines("", f)
      writeLines(paste("<r> metagenes:" ,"", round(mean(cor(t(metadata[spot.metagenes,]))), 2), sep=";"), f)

      if (length(spot.genes) < 2000)
      {
        writeLines(paste("<r> genes:" ,"", round(mean(cor(t(indata[spot.genes,]))), 2), sep=";"), f)
      }

      writeLines("", f)
      writeLines(paste("<FC> =", round(mean(indata[spot.genes,m]), 2))  ,f)
      writeLines(paste("<t-score> =", round(mean(t.g.m[spot.genes,m]), 2))  ,f)
      writeLines(paste("<p-value> =", round(10 ^ mean(log10(p.g.m[spot.genes,m])), 2))  ,f)
      writeLines(paste("<fdr> =", round(mean(fdr.g.m[spot.genes,m]), 2))  ,f)

      writeLines("", f);  writeLines("", f);  writeLines("", f)
      writeLines("Gene Statistics", f)
      writeLines("", f)
      write.csv2(out, file=f, row.names=F)

      close(f)
    }
  }

  #### Gene Set Lists ####
  if (preferences$geneset.analysis)
  {
    util.info("Writing:", file.path(dirnames["set"], "*.csv"))

    for (m in 1:ncol(indata))
    {
      gs.info <- GS.infos.samples[[m]]$GSZ.score

      pos.gs.info <- round(sort(gs.info[which(gs.info>0)],decreasing=T), 8)
      neg.gs.info <- round(sort(gs.info[which(gs.info<0)],decreasing=F), 8)

      pos.gs.info <- c(pos.gs.info, rep(0, max(0, length(neg.gs.info) - length(pos.gs.info))))
      neg.gs.info <- c(neg.gs.info, rep(0, max(0, length(pos.gs.info) - length(neg.gs.info))))

      gs.info <- data.frame(Rank=c(1:length(pos.gs.info)),
                            Upregulated=names(pos.gs.info),
                            GSZ=pos.gs.info,
                            "."=rep("",length(pos.gs.info)),
                            Downregulated=names(neg.gs.info),
                            "GSZ."=neg.gs.info)

      basename <- paste(colnames(indata)[m], ".csv", sep="")
      write.csv2(gs.info, file.path(dirnames["set"], basename), row.names=F)
    }
  }
}
