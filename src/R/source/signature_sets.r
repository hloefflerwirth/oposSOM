



  


  spotdefined.genesets = list()

  for (i in 1:length(GS.infos.overexpression$spots))
  {
    r = sapply(GS.infos.overexpression$spots[[i]]$genes, function(x)
    { 
      gene = indata[x,]
      metagene = metadata[som.nodes[x],]
  
      cor(gene, metagene)    
    })
        
    spotdefined.genesets[[names(GS.infos.overexpression$spots)[i]]] = unique(na.omit(gene.ids[GS.infos.overexpression$spots[[i]]$genes[which(r > 0.8)]]))
  }

  cat("Genes in signature sets extracted from spots:\n\n")
  print(sapply(spotdefined.genesets, length))
  

  
  
  out = cbind(Name=names(spotdefined.genesets), Type=rep(files.name ,length(spotdefined.genesets)), Literature=rep("",length(spotdefined.genesets)), Genes=sapply(spotdefined.genesets, paste, collapse=","))

  write.table(out, paste(files.name, "- Results/3rd lvl Spot Analysis/Signature Sets.csv"), sep=";", row.names=F, col.names=T)
  










