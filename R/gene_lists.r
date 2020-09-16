pipeline.geneLists <- function()
{
  
  #### Gene Localization Table ####
  o <- order(som.result$feature.BMU)
  out <- data.frame(ID=rownames(indata)[o],
                    Symbol=gene.info$names[o],
                    MeanExpression=indata.gene.mean[o],
                    Metagene=gene.info$coordinates[o],
                    Chromosome=paste( gene.info$chr.name[rownames(indata)[o]], gene.info$chr.band[rownames(indata)[o]]),
                    Description=gene.info$descriptions[o])
  
  filename <- file.path(paste(files.name, "- Results"), "CSV Sheets", "Gene localization.csv")
  util.info("Writing:", filename)
  csv.function(out, filename, row.names=FALSE)
  
	
	#### Sample GSZ Table ####
  filename <- file.path( output.paths["CSV"], "Sample GSZ scores.csv")
  util.info("Writing:", filename)
  csv.function(samples.GSZ.scores, filename)
	
	
  if(ncol(indata) < 1000)
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
  
    genes.spot.assoc <- rep("", nrow(indata) )
    names(genes.spot.assoc) <- rownames(indata)
    
    spot.list <- get(paste("spot.list.",preferences$standard.spot.modules,sep=""))
    for( i in seq_along(spot.list$spots) ) genes.spot.assoc[ spot.list$spots[[i]]$genes ] <- names(spot.list$spots)[i]
    
    for (m in 1:ncol(indata))
    {
      o <- order(p.g.m[,m])
  
      out <- data.frame(Rank=c(1:nrow(indata)),
                        ID=rownames(indata)[o],
                        Symbol=gene.info$names[o])
  
      out <- cbind(out,
                   logFC=indata[o, m],
                   WAD=WAD.g.m[o, m],
                   T.Score=t.g.m[o, m],
                   p.value=paste(p.g.m[o, m],"     ."),
                   fdr=paste(fdr.g.m[o, m],"     ."),
                   Fdr=paste(Fdr.g.m[o, m],"     ."),
                   Metagene=gene.info$coordinates[o],
                   Spot=genes.spot.assoc[o],
                   Chromosome=paste( gene.info$chr.name[rownames(indata)[o]], gene.info$chr.band[rownames(indata)[o]]),
                   Description=gene.info$descriptions[o])
  
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
      csv.function(out, file=f, row.names=FALSE)
  
      close(f)
    }
  
  
    #### Gene Set Lists ####
    if (preferences$activated.modules$geneset.analysis)
    {
      util.info("Writing:", file.path(dirnames["set"], "*.csv"))
  
      for (m in 1:ncol(indata))
      {
        gs.gsz <- spot.list.samples[[m]]$GSZ.score
  
        pos.gs.gsz <- round(sort(gs.gsz[which(gs.gsz>0)],decreasing=TRUE), 8)
        neg.gs.gsz <- round(sort(gs.gsz[which(gs.gsz<0)],decreasing=FALSE), 8)
  
        pos.gs.p <- spot.list.samples[[m]]$GSZ.p.value[names(pos.gs.gsz)]  
        neg.gs.p <- spot.list.samples[[m]]$GSZ.p.value[names(neg.gs.gsz)]
        
        pos.gs.fdr <- rep("",length(pos.gs.gsz))
        neg.gs.fdr <- rep("",length(neg.gs.gsz))

        if(ncol(indata)<100)
        {
          pos.gs.fdr <- fdrtool(pos.gs.p,statistic="pvalue",verbose=FALSE,plot=FALSE)$lfdr          
          neg.gs.fdr <- fdrtool(neg.gs.p,statistic="pvalue",verbose=FALSE,plot=FALSE)$lfdr        
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
                              "p.value"=paste(pos.gs.p,"     ."),
                              "fdr"=paste(pos.gs.fdr,"     ."),
                              "."=rep("",length(pos.gs.gsz)),
                              "Downregulated"=names(neg.gs.gsz),
                              "GSZ."=neg.gs.gsz,
                              "p.value."=paste(neg.gs.p,"     ."),
                              "fdr."=paste(neg.gs.fdr,"     ."))
  
        basename <- paste(make.names(colnames(indata)[m]), ".csv", sep="")
        csv.function(gs.info, file.path(dirnames["set"], basename), row.names=FALSE)
      }
    }
    
  }
  
}
