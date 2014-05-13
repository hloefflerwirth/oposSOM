pipeline.cancerHallmarks <- function()
{
  hallmark.sets.names <- list(
    'Angiogenesis'=c(
      "positive regulation of angiogenesis",
      "blood vessel development",
      "sprouting angiogenesis"
    ),
    'Contr. Genomic Instability'=c(
      "DNA repair",
      "response to DNA damage stimulus",
      "response to UV",
      "response to ionizing radiation"
    ),
    'Glucose Energetics'=c(
      "glycolysis",
      "positive regulation of glycolysis"
    ),
    'Inflammation'=c(
      "inflammatory response",
      "positive regulation of inflammatory response"
    ),
    'Invasion and Metastasis'=c(
      "ALONSO_METASTASIS_EMT_UP",
      "CROMER_METASTASIS_UP",
      "ODONNELL_METASTASIS_UP",
      "JAEGER_METASTASIS_UP",
      "VICENT_METASTASIS_UP",
      "BIDUS_METASTASIS_UP",
      "RAMASWAMY_METASTASIS_UP",
      "ALONSO_METASTASIS_UP",
      "negative regulation of cell adhesion",
      "negative regulation of cell-cell adhesion",
      "negative regulation of cell-matrix adhesion",
      "cell motility"
    ),
    'Proliferation'=c(
      "positive regulation of cell proliferation",
      "positive regulation of cell growth",
      "cell cycle",
      "DNA replication"
    ),
    'Replicative Immortality'=c(
      "telomere maintenance",
      "chromosome, telomeric region",
      "REACTOME_EXTENSION_OF_TELOMERES",
      "DAIRKEE_TERT_TARGETS_UP",
      "KANG_IMMORTALIZED_BY_TERT_UP"
    ),
    'Resisting Death'=c(
      "negative regulation of apoptotic process",
      "negative regulation of autophagy"
    )
  )

  hallmark.sets.genes <- lapply(hallmark.sets.names, function(x)
  {
    unique(unlist(sapply(gs.def.list[x], function(y) { y$Genes })))
  })

  hallmark.sets.ids <- lapply(hallmark.sets.genes, function(x)
  {
    unique(names(gene.ids)[which(gene.ids %in% x)])
  })

  ### GSZ profiles

  hallmark.sets.list <- lapply(hallmark.sets.genes, function(x) list(Genes=x,Type=""))

  cl <- makeCluster(preferences$max.parallel.cores)

  hallmark.GSZ.matrix <- sapply(1:ncol(indata), function(m)
  {
    return(GeneSet.GSZ(unique.protein.ids, batch.t.g.m[, m], hallmark.sets.list, sort=F, cluster=cl))
  })

  hallmark.spot.enrichment <- sapply(spot.list.overexpression$spots, function(x)
  {
    spot.ens.ids <- unique(na.omit(gene.ids[x$genes]))
    return(GeneSet.Fisher(spot.ens.ids, unique.protein.ids, hallmark.sets.list, sort=F, cluster=cl))
  })

  try({ stopCluster(cl) }, silent=T)

  hallmark.sets.profiles <- list()

  for (i in 1:nrow(hallmark.GSZ.matrix))
  {
    hallmark.sets.profiles[[rownames(hallmark.GSZ.matrix)[i]]] <- hallmark.GSZ.matrix[i,]
    names(hallmark.sets.profiles[[rownames(hallmark.GSZ.matrix)[i]]]) <- colnames(indata)
  }

  ### Output
  filename <- file.path(paste(files.name, "- Results"), "Geneset Analysis", "0verview Cancer Hallmarks.pdf")
  util.info("Writing:", filename)
  pdf(filename, 21/2.54, 29.7/2.54)

  layout(matrix(c(1:8), 4, byrow=T), widths=c(3, 1))

  for (i in 1:length(hallmark.sets.profiles))
  {
    hallmark.sets.group.profiles <-
      tapply(hallmark.sets.profiles[[i]], group.labels, c)[unique(group.labels)]

    par(mar=c(5,3,4,2))

    boxplot(hallmark.sets.group.profiles, col=unique.group.colors, las=2,
            main=names(hallmark.sets.names)[i], ylim=range(hallmark.GSZ.matrix))

    abline(h=0, lty=2)

    n.map <- matrix(0,preferences$dim.1stLvlSom,preferences$dim.1stLvlSom)
    gs.nodes <- som.nodes[hallmark.sets.ids[[i]]]
    n.map[as.numeric(names(table(gs.nodes)))] <- table(gs.nodes)
    n.map[which(n.map==0)] <- NA
    n.map <- matrix(n.map, preferences$dim.1stLvlSom)

    par(mar=c(5,1,4,1))

    lim <- c(1,preferences$dim.1stLvlSom) + preferences$dim.1stLvlSom*0.01*c(-1,1)
    colr <- colramp(1000)[(na.omit(as.vector(n.map)) - min(n.map,na.rm=T)) /
                          max(1, (max(n.map,na.rm=T) - min(n.map,na.rm=T))) *
                          999 + 1]

    plot(which(!is.na(n.map), arr.ind=T), xlim=lim, ylim=lim, pch=16,
         axes=F, xlab="",ylab="", xaxs="i", yaxs="i", col=colr,
         cex=0.5 + na.omit(as.vector(n.map)) / max(n.map,na.rm=T) * 2.8)

    title(sub=paste("# features =", length(hallmark.sets.ids[[i]]), ", max =",
                    max(n.map,na.rm=T)),line=0)

    box()
  }

  par(mfrow=c(1,1))

  hallmark.spot.enrichment.bin <- hallmark.spot.enrichment
  hallmark.spot.enrichment.bin[which(hallmark.spot.enrichment >= 0)] <- 1
  hallmark.spot.enrichment.bin[which(hallmark.spot.enrichment < 0.1)] <- 0.75
  hallmark.spot.enrichment.bin[which(hallmark.spot.enrichment < 0.01)] <- 0.5
  hallmark.spot.enrichment.bin[which(hallmark.spot.enrichment < 0.001)] <- 0.25
  hallmark.spot.enrichment.bin[which(hallmark.spot.enrichment < 0.00001)] <- 0.0

  heatmap(hallmark.spot.enrichment.bin,
          col=c("black","red3","orange","yellow","gray92"),
          Colv=NA, Rowv=NA, scale="n", mar=c(10,10),
          main="Hallmark spot enrichment", margins=c(5,13))

  par(new=T, mar=c(5,0,0,3))
  frame()

  legend("bottomright", c("p < 0.1","p < 0.01","p < 0.001","p < 0.00001"),
         pch=15, col=rev(c("black","red3","orange","yellow")))

  dev.off()
}
