pipeline.geneLists <- function(env)
{
  
  #### Gene Localization Table ####
  
  spot.list <- env[[paste("spot.list.", env$preferences$standard.spot.modules,sep="")]]
  gene.module.info <- rep("",nrow(env$indata))
  names(gene.module.info) <- rownames(env$indata)
  for( x in names(spot.list$spots) )
  {
    gene.module.info[ spot.list$spots[[x]]$genes ] <- x
  }
  
  o <- order(env$som.result$feature.BMU)
  out <- data.frame(ID=rownames(env$indata)[o],
                    Symbol=env$gene.info$names[o],
                    MeanExpression=env$indata.gene.mean[o],
                    Metagene=env$gene.info$coordinates[o],
                    Module=gene.module.info[o],
                    Chromosome=paste( env$gene.info$chr.name[rownames(env$indata)[o]], env$gene.info$chr.band[rownames(env$indata)[o]]),
                    Description=env$gene.info$descriptions[o])
  
  filename <- file.path("CSV Sheets", "Gene localization.csv")
  util.info("Writing:", filename)
  env$csv.function(out, filename, row.names=FALSE)
  
	
	#### Sample GSZ Table ####
  filename <- file.path( env$output.paths["CSV"], "Sample GSZ scores.csv")
  util.info("Writing:", filename)
  env$csv.function(env$samples.GSZ.scores, filename)
	
	
  #### sample-wise gene and gene set tables ####
  
  # if(ncol(env$indata) < 1000)
  # {
  #   dirnames <- c("global"=file.path(env$output.paths["CSV"], "Gene Lists - Global"),
  # #                "local"=file.path(output.paths["CSV"], "Gene Lists - Local"),
  #                 "set"=file.path(env$output.paths["CSV"], "Gene Set Lists - Global"))
  # 
  #   for (dirname in dirnames)
  #   {
  #     dir.create(dirname, showWarnings=FALSE)
  #   }
  # 
  #   #### Global Gene Lists ####
  #   util.info("Writing:", file.path(dirnames["global"], "*.csv"))
  # 
  #   genes.spot.assoc <- rep("", nrow(env$indata) )
  #   names(genes.spot.assoc) <- rownames(env$indata)
  #   
  #   spot.list <- env[[paste("spot.list.",env$preferences$standard.spot.modules,sep="")]]
  #   for( i in seq_along(spot.list$spots) ) genes.spot.assoc[ spot.list$spots[[i]]$genes ] <- names(spot.list$spots)[i]
  #   
  #   for (m in 1:ncol(env$indata))
  #   {
  #     o <- order(env$p.g.m[,m])
  # 
  #     out <- data.frame(Rank=c(1:nrow(env$indata)),
  #                       ID=rownames(env$indata)[o],
  #                       Symbol=env$gene.info$names[o])
  # 
  #     out <- cbind(out,
  #                  logFC=env$indata[o, m],
  #                  p.value=env$p.g.m[o, m],
  #                  fdr=env$fdr.g.m[o, m],
  #                  Metagene=env$gene.info$coordinates[o],
  #                  Spot=genes.spot.assoc[o],
  #                  Chromosome=paste( env$gene.info$chr.name[rownames(env$indata)[o]], env$gene.info$chr.band[rownames(env$indata)[o]]),
  #                  Description=env$gene.info$descriptions[o])
  # 
  #     basename <- paste(make.names(colnames(env$indata)[m]), ".csv", sep="")
  #     f <- file(file.path(dirnames["global"], basename), "w")
  # 
  #     writeLines("Sample Summary:", f)
  #     writeLines("", f)
  #     writeLines(paste("%DE:" ,"", round(env$perc.DE.m[colnames(env$indata)[m]], 2), sep=";"), f)
  #     writeLines(paste("#genes with fdr < 0.2" ,"", length(which(env$fdr.g.m[,m] < 0.2)) , sep=";"), f)
  #     writeLines(paste("#genes with fdr < 0.1" ,"", length(which(env$fdr.g.m[,m] < 0.1)) , sep=";"), f)
  #     writeLines(paste("#genes with fdr < 0.05" ,"", length(which(env$fdr.g.m[,m] < 0.05)) , sep=";"), f)
  #     writeLines(paste("#genes with fdr < 0.01" ,"", length(which(env$fdr.g.m[,m] < 0.01)) , sep=";"), f)
  #     writeLines("", f)
  #     writeLines(paste("<FC> =", round(mean(env$indata[,m]), 2))  ,f)
  #     writeLines(paste("<p-value> =", round(10 ^ mean(log10(env$p.g.m[,m])), 2))  ,f)
  #     writeLines(paste("<fdr> =", round(mean(env$fdr.g.m[,m]), 2))  ,f)
  #     writeLines("", f);  writeLines("", f);  writeLines("", f)
  #     writeLines("Gene Statistics", f)
  #     writeLines("", f)
  #     env$csv.function(out, file=f, row.names=FALSE)
  # 
  #     close(f)
  #   }
  # 
  # 
  #   #### Gene Set Lists ####
  #   if (env$preferences$activated.modules$geneset.analysis)
  #   {
  #     util.info("Writing:", file.path(dirnames["set"], "*.csv"))
  # 
  #     for (m in 1:ncol(env$indata))
  #     {
  #       gs.gsz <- env$spot.list.samples[[m]]$GSZ.score
  # 
  #       pos.gs.gsz <- round(sort(gs.gsz[which(gs.gsz>0)],decreasing=TRUE), 8)
  #       neg.gs.gsz <- round(sort(gs.gsz[which(gs.gsz<0)],decreasing=FALSE), 8)
  # 
  #       pos.gs.p <- env$spot.list.samples[[m]]$GSZ.p.value[names(pos.gs.gsz)]  
  #       neg.gs.p <- env$spot.list.samples[[m]]$GSZ.p.value[names(neg.gs.gsz)]
  #       
  #       pos.gs.fdr <- rep("",length(pos.gs.gsz))
  #       neg.gs.fdr <- rep("",length(neg.gs.gsz))
  # 
  #       if(ncol(env$indata)<100)
  #       {
  #         pos.gs.fdr <- fdrtool(pos.gs.p,statistic="pvalue",verbose=FALSE,plot=FALSE)$lfdr          
  #         neg.gs.fdr <- fdrtool(neg.gs.p,statistic="pvalue",verbose=FALSE,plot=FALSE)$lfdr        
  #       }
  #             
  #       pos.gs.gsz <- c(pos.gs.gsz, rep(0, max(0, length(neg.gs.gsz) - length(pos.gs.gsz))))
  #       neg.gs.gsz <- c(neg.gs.gsz, rep(0, max(0, length(pos.gs.gsz) - length(neg.gs.gsz))))
  #       pos.gs.p <- c(pos.gs.p, rep(1, max(0, length(neg.gs.p) - length(pos.gs.p))))
  #       neg.gs.p <- c(neg.gs.p, rep(1, max(0, length(pos.gs.p) - length(neg.gs.p))))
  #       pos.gs.fdr <- c(pos.gs.fdr, rep(1, max(0, length(neg.gs.fdr) - length(pos.gs.fdr))))
  #       neg.gs.fdr <- c(neg.gs.fdr, rep(1, max(0, length(pos.gs.fdr) - length(neg.gs.fdr))))    
  #       
  #       gs.info <- data.frame("Rank"=c(seq_along(pos.gs.gsz)),
  #                             "Upregulated"=names(pos.gs.gsz),
  #                             "GSZ"=pos.gs.gsz,
  #                             "p.value"=paste(pos.gs.p,"     ."),
  #                             "fdr"=paste(pos.gs.fdr,"     ."),
  #                             "."=rep("",length(pos.gs.gsz)),
  #                             "Downregulated"=names(neg.gs.gsz),
  #                             "GSZ."=neg.gs.gsz,
  #                             "p.value."=paste(neg.gs.p,"     ."),
  #                             "fdr."=paste(neg.gs.fdr,"     ."))
  # 
  #       basename <- paste(make.names(colnames(env$indata)[m]), ".csv", sep="")
  #       env$csv.function(gs.info, file.path(dirnames["set"], basename), row.names=FALSE)
  #     }
  #   }
  #   
  # }
  
}
