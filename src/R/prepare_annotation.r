pipeline.prepareAnnotation <- function()
{
  gene.ids <<- rep("", nrow(indata))
  names(gene.ids) <<- rownames(indata)

  gene.names <<- rep("", nrow(indata))
  names(gene.names) <<- rownames(indata)

  gene.descriptions <<- rep("", nrow(indata))
  names(gene.descriptions) <<- rownames(indata)

  gene.positions <<- rep("", nrow(indata))
  names(gene.positions) <<- rownames(indata)
  gene.positions.table <<- matrix(0,0,0)
  gene.positions.list <<- list()

  if (preferences$ensembl.dataset == "" || preferences$ensembl.rowname.ids == "")
  {
    preferences$geneset.analysis <<- F
    return()
  }

  require.bioconductor("biomaRt")

  mart<-useMart('ensembl')
  mart<-useDataset(preferences$ensembl.dataset, mart=mart)

  biomart.table = getBM(c(preferences$ensembl.rowname.ids,
                          "external_gene_id",
                          "description",
                          "ensembl_gene_id",
                          "chromosome_name",
                          "band"),
                        preferences$ensembl.rowname.ids,
                        rownames(indata), mart, checkFilters=F)


  if (nrow(biomart.table) == 0)
  {
    preferences$geneset.analysis <<- F
    util.warn("Could not resolve rownames. Possibly wrong ensembl.rowname.ids")
  } else
  {
    h <- biomart.table[,2]
    names(h) <- biomart.table[,1]
    gene.names[as.character(unique(biomart.table[,1]))] <<- h[as.character(unique(biomart.table[,1]))]

    h <- biomart.table[,3]
    names(h) <- biomart.table[,1]
    gene.descriptions[as.character(unique(biomart.table[,1]))] <<- h[as.character(unique(biomart.table[,1]))]

    h <- biomart.table[,4]
    names(h) <- biomart.table[,1]
    gene.ids[as.character(unique(biomart.table[,1]))] <<- h[as.character(unique(biomart.table[,1]))]
    gene.ids <<- gene.ids[which(gene.ids != "")]
    gene.ids <<- gene.ids[which(names(gene.ids) %in% rownames(indata))]

    h <- paste(biomart.table[,"chromosome_name"], gsub("\\..*$","", biomart.table[,"band"]))
    names(h) <- biomart.table[,1]
    gene.positions[as.character(unique(biomart.table[,1]))] <<- h[as.character(unique(biomart.table[,1]))]
    gene.positions <<- gene.positions[which(gene.positions != "")]

    gene.positions.table <<- do.call(rbind, strsplit(gene.positions, " "))

    junk.chrnames <- names(which(table(gene.positions.table[,1]) < 20))
    junk.chrnames <- union(junk.chrnames, names(which(tapply(gene.positions.table[,2], gene.positions.table[,1], function(x) { length(unique(x)) }) == 1)))
    gene.positions <<- gene.positions[which(!gene.positions.table[,1] %in% junk.chrnames)]
    gene.positions.table <<- gene.positions.table[which(!gene.positions.table[,1] %in% junk.chrnames),]

    gene.positions.list <<- tapply(rownames(gene.positions.table), gene.positions.table[,1], c)
    gene.positions.list <<- lapply(gene.positions.list, function(x) { tapply(x, gene.positions.table[x,2], c) })
  }

  if (!preferences$geneset.analysis) {
    return()
  }

  unique.protein.ids <<- unique(gene.ids)

  biomart.table <- getBM(c("ensembl_gene_id", "go_id", "name_1006", "namespace_1003"), "ensembl_gene_id", unique.protein.ids, mart, checkFilters=F)
  gs.def.list <<- tapply(biomart.table[,1], biomart.table[,2], c)
  gs.def.list <<- lapply(gs.def.list, function(x) { list(Genes=x, Type="") })
  gs.def.list <<- gs.def.list[- which(names(gs.def.list) == "")]


  if (length(gs.def.list) > 0)
  {
    util.info("Download of", length(gs.def.list), "GO sets with", sum(sapply(sapply(gs.def.list, head, 1), length)), "entries")

    ## simple small-gs-filtering
    gs.def.list <<- gs.def.list[- which(sapply(sapply(gs.def.list, head, 1), length) < 10)]
    util.info("Filtered to", length(gs.def.list), "GO sets with", sum(sapply(sapply(gs.def.list, head, 1), length)), "entries")

    biomart.table[,4] <- sub("biological_process", "BP", biomart.table[,4])
    biomart.table[,4] <- sub("molecular_function", "MF", biomart.table[,4])
    biomart.table[,4] <- sub("cellular_component", "CC", biomart.table[,4])

    for (i in 1:length(gs.def.list))
    {
      o <- match(names(gs.def.list)[i], biomart.table[,2])
      names(gs.def.list)[i] <<- biomart.table[o, 3]
      gs.def.list[[i]]$Type <<- biomart.table[o, 4]
    }

    gs.def.list <<- gs.def.list[order(names(gs.def.list))]
  } else
  {
    util.warn("No GO annotation found")
  }

  chr.gs.list <- lapply(gene.positions.list, function(x) { list(Genes=gene.ids[unlist(x)], Type="Chr") })
  names(chr.gs.list) <- paste("Chr", names(gene.positions.list))
  gs.def.list <<- c(gs.def.list, chr.gs.list)

  small.gs <- which(sapply(sapply(gs.def.list, head, 1), length) < 10)

  if (length(small.gs) > 0)
  {
    gs.def.list <<- gs.def.list[- small.gs]
  }

  if (length(preferences$geneset.custom.list) > 0)
  {
    gs.def.list <<- c(gs.def.list, preferences$geneset.custom.list)
  }

  gs.def.list <<- lapply(gs.def.list, function(x) {
    x$Genes <- intersect(x$Genes, unique.protein.ids)
    return(x)
  })

  small.gs <- which(sapply(sapply(gs.def.list, head, 1), length) < 2)

  if (length(small.gs) > 0)
  {
    gs.def.list <<- gs.def.list[- small.gs]
  }

  if (length(gs.def.list) > 0)
  {
    gs.def.list <<- lapply(gs.def.list, function(x) { names(x$Genes) = NULL; return(x) })
    gs.def.list.categories <<- sapply(gs.def.list, function(x) { x$Type })
    util.info("In total", length(gs.def.list), "gene sets to be considered in analysis")
  } else
  {
    preferences$geneset.analysis <<- F
    util.warn("No Geneset information -> turning off GS analysis")
  }
}
