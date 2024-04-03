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
  filename <- file.path( "CSV Sheets", "Sample GSZ scores.csv")
  util.info("Writing:", filename)
  env$csv.function(env$samples.GSZ.scores, filename)
  
}
