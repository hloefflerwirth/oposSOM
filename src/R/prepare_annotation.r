pipeline.prepareAnnotation <- function()
{
  empty.vec.chr <- rep("", nrow(indata))
  names(empty.vec.chr) <- rownames(indata)
  empty.vec.num <- rep(NA, nrow(indata))
  names(empty.vec.num) <- rownames(indata)
  mode(empty.vec.num)  <- "numeric"
  
  gene.info <<- list( ids = empty.vec.chr, 
                      names = empty.vec.chr,
                      descriptions = empty.vec.chr,
                      chr.name = empty.vec.chr,
                      chr.band = empty.vec.chr,
                      chr.start = empty.vec.num )
  
  chromosome.list <<- list()

  if (!util.call(biomart.available, environment()))
  {
    util.warn("biomaRt seems to be down.")
    util.warn("Disabling geneset analysis.")
    preferences$activated.modules$geneset.analysis <<- FALSE
    return()
  }

  if (!preferences$database.dataset %in% c("auto", ""))
  {
    biomart.table <- NULL

    try({
      mart <- useMart(biomart=preferences$database.biomart, host=preferences$database.host)
      mart <- useDataset(preferences$database.dataset, mart=mart)

      query = c("hgnc_symbol","wikigene_name","uniprot_genename")[ which( c("hgnc_symbol","wikigene_name","uniprot_genename") %in% listAttributes(mart)[,1] ) ][1]
      suppressWarnings({  biomart.table <-
        getBM(c(preferences$database.id.type, query),
              preferences$database.id.type,
              rownames(indata)[seq(1,nrow(indata),length.out=100)],
              mart, checkFilters=FALSE)  })
    }, silent=TRUE)

    if (is.null(biomart.table) || nrow(biomart.table) == 0)
    {
      util.warn("Invalid annotation parameters. Trying autodetection...")
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
    preferences$activated.modules$geneset.analysis <<- FALSE
    return()
  }
  

  mart <- useMart(biomart=preferences$database.biomart, host=preferences$database.host)
  mart <- useDataset(preferences$database.dataset, mart=mart)

  query = c("wikigene_name","hgnc_symbol","uniprot_genename")[ which( c("wikigene_name","hgnc_symbol","uniprot_genename") %in% listAttributes(mart)[,1] ) ][1]
  suppressWarnings({  biomart.table <- getBM(c(preferences$database.id.type,
                           query,"description","ensembl_gene_id",
                           "chromosome_name","band","start_position"),
                         preferences$database.id.type,
                         rownames(indata), mart, checkFilters=FALSE)  })

  biomart.table <- biomart.table[ which(biomart.table[,1]%in%rownames(indata)), ]
  
  if (nrow(biomart.table) == 0)
  {
    util.warn("Could not resolve rownames. Possibly wrong database.id.type")
    util.warn("Disabling geneset analysis.")
    preferences$activated.modules$geneset.analysis <<- FALSE
  } else
  {
    h <- biomart.table[,2]
    names(h) <- biomart.table[,1]
    gene.info$names[as.character(unique(biomart.table[,1]))] <<- h[as.character(unique(biomart.table[,1]))]

    h <- biomart.table[,"description"]
    names(h) <- biomart.table[,1]
    gene.info$descriptions[as.character(unique(biomart.table[,1]))] <<- h[as.character(unique(biomart.table[,1]))]

    h <- biomart.table[,"ensembl_gene_id"]
    names(h) <- biomart.table[,1]
    gene.info$ids[as.character(unique(biomart.table[,1]))] <<- h[as.character(unique(biomart.table[,1]))]

    h <- biomart.table[,"chromosome_name"]
    names(h) <- biomart.table[,1]
    gene.info$chr.name[as.character(unique(biomart.table[,1]))] <<- h[as.character(unique(biomart.table[,1]))]
    
    h <- gsub("\\..*$","", biomart.table[,"band"])
    names(h) <- biomart.table[,1]
    gene.info$chr.band[as.character(unique(biomart.table[,1]))] <<- h[as.character(unique(biomart.table[,1]))]
    
    h <- biomart.table[,"start_position"]
    names(h) <- biomart.table[,1]
    gene.info$chr.start[as.character(unique(biomart.table[,1]))] <<- h[as.character(unique(biomart.table[,1]))]
    

    gene.positions.table <- cbind( gene.info$chr.name, gene.info$chr.band )
    gene.positions.table <- gene.positions.table[ which( gene.positions.table[,1] != "" & gene.positions.table[,2] != "" ) , ]
    skip.chrnames <- names(which(table(gene.info$chr.name) < 20)) # filter low abundand chromosome information
    gene.positions.table <- gene.positions.table[ which( !gene.positions.table[,1] %in% skip.chrnames ) , ]
      
    if(nrow(gene.positions.table)>0)
    {
      chromosome.list <<- tapply(rownames(gene.positions.table), gene.positions.table[,1], c)
      chromosome.list <<- lapply(chromosome.list, function(x) { tapply(x, gene.positions.table[x,2], c) })
    }

  }

  if (!preferences$activated.modules$geneset.analysis) {
    return()
  }

  unique.protein.ids <<- unique(gene.info$ids)
  unique.protein.ids <<- unique.protein.ids[which(unique.protein.ids!="")]
  
  suppressWarnings({  biomart.table <- getBM(c("ensembl_gene_id", "go_id", "name_1006", "namespace_1003"), "ensembl_gene_id", unique.protein.ids, mart, checkFilters=FALSE) })
  biomart.table <- biomart.table[which( apply(biomart.table,1,function(x) sum(x=="") ) == 0 ),]  
  gs.def.list <<- tapply(biomart.table[,1], biomart.table[,2], c)
  gs.def.list <<- lapply(gs.def.list, function(x) { list(Genes=x, Type="") })


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

  if(length(chromosome.list)>0)
  {
    chr.gs.list <- lapply(chromosome.list, function(x) { list(Genes=gene.info$ids[unlist(x)], Type="Chr") })
    names(chr.gs.list) <- paste("Chr", names(chromosome.list))
    gs.def.list <<- c(gs.def.list, chr.gs.list)
  }


  # load custom genesets
  data(opossom.genesets)
  gs.def.list  <<- c(gs.def.list, opossom.genesets)

  gs.def.list <<- lapply(gs.def.list, function(x) {
    x$Genes <- intersect(x$Genes, unique.protein.ids)
    return(x)
  })


  gs.def.list <<- gs.def.list[ which(sapply(sapply(gs.def.list, head, 1), length) >= 2) ]

  if (length(gs.def.list) > 0)
  {
    gs.def.list <<- lapply(gs.def.list, function(x) { names(x$Genes) = NULL; return(x) })
    util.info("In total", length(gs.def.list), "gene sets to be considered in analysis")
  } else
  {
    preferences$activated.modules$geneset.analysis <<- FALSE
    util.warn("No Geneset information -> turning off GS analysis")
  }
}
