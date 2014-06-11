pipeline.signatureSets <- function()
{
  spotdefined.genesets <- list()

  for (i in 1:length(spot.list.overexpression$spots))
  {
    r <- sapply(spot.list.overexpression$spots[[i]]$genes, function(x)
    {
      gene <- indata[x,]
      metagene <- metadata[som.nodes[x],]
      cor(gene, metagene)
    })

    spotdefined.genesets[[names(spot.list.overexpression$spots)[i]]] <-
      unique(na.omit(gene.ids[spot.list.overexpression$spots[[i]]$genes[which(r > 0.8)]]))
  }

  out <- cbind(Name=names(spotdefined.genesets),
               Type=rep(files.name ,length(spotdefined.genesets)),
               Literature=rep("",length(spotdefined.genesets)),
               Genes=sapply(spotdefined.genesets, paste, collapse=","))

  filename <- file.path(paste(files.name, "- Results"),
                        "3rd lvl Spot Analysis",
                        "Signature Sets.csv")

  util.info("Writing:", filename)
  write.table(out, filename, sep=";", row.names=FALSE, col.names=TRUE)
}
