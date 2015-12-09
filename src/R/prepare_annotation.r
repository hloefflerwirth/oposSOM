pipeline.prepareAnnotation <- function()
{
  util.info("Loading gene annotations from database. This may take several minutes until next notification.")
  
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

  if (!biomart.available())
  {
    util.warn("biomaRt seems to be down.")
    util.warn("Disabling geneset analysis.")
    preferences$geneset.analysis <<- FALSE
    return()
  }

  if (!preferences$database.dataset %in% c("auto", ""))
  {
    biomart.table <- NULL

    try({
      mart <- useMart('ENSEMBL_MART_ENSEMBL',host="www.ensembl.org")
      mart <- useDataset(preferences$database.dataset, mart=mart)

      query = c("hgnc_symbol","wikigene_name","uniprot_genename")[ which( c("hgnc_symbol","wikigene_name","uniprot_genename") %in% listAttributes(mart)[,1] ) ][1]
      biomart.table <-
        getBM(c(preferences$database.id.type, query),
              preferences$database.id.type,
              rownames(indata)[seq(1,nrow(indata),length.out=100)],
              mart, checkFilters=FALSE)
    }, silent=TRUE)

    if (is.null(biomart.table) || nrow(biomart.table) == 0)
    {
      util.warn("Invalid annotation parameters. Try to autodetect...")
      preferences$database.dataset <<- "auto"
    }
  }

  if (preferences$database.dataset == "auto")
  {
    util.call(pipeline.detectEnsemblDataset, environment())
  }

  if (preferences$database.dataset == "" || preferences$database.id.type == "")
  {
    util.warn("Could not find valid annotation parameters.")
    util.warn("Disabling geneset analysis.")
    preferences$geneset.analysis <<- FALSE
    return()
  }
  

  mart <- useMart('ENSEMBL_MART_ENSEMBL',host="www.ensembl.org")
  mart <- useDataset(preferences$database.dataset, mart=mart)

  query = c("wikigene_name","hgnc_symbol","uniprot_genename")[ which( c("wikigene_name","hgnc_symbol","uniprot_genename") %in% listAttributes(mart)[,1] ) ][1]
  biomart.table <- getBM(c(preferences$database.id.type,
                           query,
                           "description",
                           "ensembl_gene_id",
                           "chromosome_name",
                           "band"),
                         preferences$database.id.type,
                         rownames(indata), mart, checkFilters=FALSE)

  if (nrow(biomart.table) == 0)
  {
    util.warn("Could not resolve rownames. Possibly wrong database.id.type")
    util.warn("Disabling geneset analysis.")
    preferences$geneset.analysis <<- FALSE
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

    if(length(gene.positions)>0)
    {
      gene.positions.list <<- tapply(rownames(gene.positions.table), gene.positions.table[,1], c)
      gene.positions.list <<- lapply(gene.positions.list, function(x) { tapply(x, gene.positions.table[x,2], c) })
    }
  }

  if (!preferences$geneset.analysis) {
    return()
  }

  unique.protein.ids <<- unique(gene.ids)

  biomart.table <- getBM(c("ensembl_gene_id", "go_id", "name_1006", "namespace_1003"), "ensembl_gene_id", unique.protein.ids, mart, checkFilters=FALSE)
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

    for (i in seq_along(gs.def.list))
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

  if(length(gene.positions.list)>0)
  {
    chr.gs.list <- lapply(gene.positions.list, function(x) { list(Genes=gene.ids[unlist(x)], Type="Chr") })
    names(chr.gs.list) <- paste("Chr", names(gene.positions.list))
    gs.def.list <<- c(gs.def.list, chr.gs.list)
  }


  # load custom genesets
  data(opossom.genesets)
  gs.def.list  <<- c(gs.def.list, opossom.genesets)

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
    sapply(gs.def.list, function(x) { x$Type }) <<- sapply(gs.def.list, function(x) { x$Type })
    util.info("In total", length(gs.def.list), "gene sets to be considered in analysis")
  } else
  {
    preferences$geneset.analysis <<- FALSE
    util.warn("No Geneset information -> turning off GS analysis")
  }
}
