

if (length(preferences$differences.list) > 0)
{
  
  for (d in 1:length(preferences$differences.list))
  {    
    samples.indata = list(preferences$differences.list[[d]][[1]], preferences$differences.list[[d]][[2]])
    samples.indata.original = list(which(colnames(indata.original) %in% colnames(indata)[samples.indata[[1]]])  ,  which(colnames(indata.original) %in% colnames(indata)[samples.indata[[2]]]))    
      
    indata = cbind(indata,   rowMeans(indata.original[,samples.indata.original[[2]],drop=F]) - rowMeans(indata.original[,samples.indata.original[[1]],drop=F]))  
    colnames(indata)[ncol(indata)] = names(preferences$differences.list)[d]
  }  
  
}
colnames(indata) = paste("logFC ", colnames(indata))



cat("|")
for (i in 1:50) cat(" ")
cat("|\n|")
flush.console()


#i=length(gs.def.list)
for (i in 1:length(gs.def.list))
{
  genes = names(gene.ids)[which(gene.ids %in% gs.def.list[[i]]$Genes)]
  
  out = data.frame(AffyID = names(gene.ids[genes]),
                    EnsemblID = gene.ids[genes],
                    Metagene = genes.coordinates[genes],
                    Max.expression.sample = colnames(indata)[apply(indata[genes, ,drop=F], 1, which.max)],
                    GeneSymbol = gene.names[genes],
                    Description = gene.descriptions[genes])

  
  out = cbind(out, indata[genes, ,drop=F])    
    

    
#  write.csv2(out, "a.csv", row.names=F)
  write.csv2(out, paste(files.name, " - Results/Geneset Analysis/", substring(make.names(names(gs.def.list)[i]), 1, 100),".csv", sep=""), row.names=F)

  
  out.intervals = round(seq(1, length(gs.def.list), length.out=50+1))[-1]
  cat(paste(rep("#",length(which(out.intervals == i))), collapse=""));  flush.console()
}

cat("|\n\n");  flush.console()