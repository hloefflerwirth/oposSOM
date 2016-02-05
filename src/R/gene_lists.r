pipeline.geneLists <- function()
{
  dirnames <- c("global"=file.path(output.paths["CSV"], "Gene Lists - Global"),
#                "local"=file.path(output.paths["CSV"], "Gene Lists - Local"),
                "set"=file.path(output.paths["CSV"], "Gene Set Lists - Global"))

  for (dirname in dirnames)
  {
    dir.create(dirname, showWarnings=FALSE)
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
                 Metagene=gene.coordinates[o],
                 Chromosome=gene.positions[rownames(indata)[o]],
                 Description=gene.descriptions[o])

    basename <- paste(make.names(colnames(indata)[m]), ".csv", sep="")
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
    write.csv2(out, file=f, row.names=FALSE)

    close(f)
  }


  #### Gene Set Lists ####
  if (preferences$geneset.analysis)
  {
    util.info("Writing:", file.path(dirnames["set"], "*.csv"))

    for (m in 1:ncol(indata))
    {
      gs.gsz <- spot.list.samples[[m]]$GSZ.score

      pos.gs.gsz <- round(sort(gs.gsz[which(gs.gsz>0)],decreasing=TRUE), 8)
      neg.gs.gsz <- round(sort(gs.gsz[which(gs.gsz<0)],decreasing=FALSE), 8)

      pos.gs.p <- rep("",length(pos.gs.gsz))
      neg.gs.p <- rep("",length(neg.gs.gsz))
      
      pos.gs.fdr <- rep("",length(pos.gs.gsz))
      neg.gs.fdr <- rep("",length(neg.gs.gsz))
      
      if(preferences$geneset.analysis.exact)
      {
        pos.gs.p <- spot.list.samples[[m]]$GSZ.p.value[names(pos.gs.gsz)]  
        neg.gs.p <- spot.list.samples[[m]]$GSZ.p.value[names(neg.gs.gsz)]
        
        if(ncol(indata)<100)
        {
          pos.gs.fdr <- fdrtool(pos.gs.p,statistic="pvalue",verbose=F,plot=F)$lfdr          
          neg.gs.fdr <- fdrtool(neg.gs.p,statistic="pvalue",verbose=F,plot=F)$lfdr        
        }
      }
            
      pos.gs.gsz <- c(pos.gs.gsz, rep(0, max(0, length(neg.gs.gsz) - length(pos.gs.gsz))))
      neg.gs.gsz <- c(neg.gs.gsz, rep(0, max(0, length(pos.gs.gsz) - length(neg.gs.gsz))))
      pos.gs.p <- c(pos.gs.p, rep(1, max(0, length(neg.gs.p) - length(pos.gs.p))))
      neg.gs.p <- c(neg.gs.p, rep(1, max(0, length(pos.gs.p) - length(neg.gs.p))))
      pos.gs.fdr <- c(pos.gs.fdr, rep(1, max(0, length(neg.gs.fdr) - length(pos.gs.fdr))))
      neg.gs.fdr <- c(neg.gs.fdr, rep(1, max(0, length(pos.gs.fdr) - length(neg.gs.fdr))))    
      
      gs.info <- data.frame("Rank"=c(seq_along(pos.gs.gsz)),
                            "Upregulated"=names(pos.gs.gsz),
                            "GSZ"=pos.gs.gsz,
                            "p.value"=pos.gs.p,
                            "fdr"=pos.gs.fdr,
                            "."=rep("",length(pos.gs.gsz)),
                            "Downregulated"=names(neg.gs.gsz),
                            "GSZ."=neg.gs.gsz,
                            "p.value."=neg.gs.p,
                            "fdr."=neg.gs.fdr)

      basename <- paste(make.names(colnames(indata)[m]), ".csv", sep="")
      write.csv2(gs.info, file.path(dirnames["set"], basename), row.names=FALSE)
    }
  }
}
