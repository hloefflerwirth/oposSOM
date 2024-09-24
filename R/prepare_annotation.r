pipeline.prepareAnnotation <- function(env)
{
  empty.vec.chr <- rep("", nrow(env$indata))
  names(empty.vec.chr) <- rownames(env$indata)
  empty.vec.num <- rep(NA, nrow(env$indata))
  names(empty.vec.num) <- rownames(env$indata)
  mode(empty.vec.num)  <- "numeric"
  
  env$gene.info <- list( ids = empty.vec.chr,
                      names = empty.vec.chr,
                      descriptions = empty.vec.chr,
                      chr.name = empty.vec.chr,
                      chr.band = empty.vec.chr,
                      chr.start = empty.vec.num,
                      ensembl.mapping = data.frame( indata.id=rownames(env$indata), ensembl.id=rep("", nrow(env$indata)), stringsAsFactors=FALSE )   )
  
  if(!env$preferences$activated.modules$primary.analysis)
  {
    env$gene.info$coordinates <- apply( env$som.result$node.summary[env$som.result$feature.BMU,c("x","y")], 1, paste, collapse=" x " )
    names(env$gene.info$coordinates) <- rownames(env$indata)
  }
  
  env$chromosome.list <- list()

  if (!biomart.available(env))
  {
    util.warn("Requested biomaRt host seems to be down.")
    util.warn("Disabling geneset analysis.")
    env$preferences$activated.modules$geneset.analysis <- FALSE
    env$preferences$activated.modules$psf.analysis <- FALSE
    return(env)
  }

  if (!env$preferences$database.dataset %in% c("auto", ""))
  {
    biomart.table <- NULL

    try({
      mart <- useMart(biomart=env$preferences$database.biomart, host=env$preferences$database.host)
      mart <- useDataset(env$preferences$database.dataset, mart=mart)

      query = c("hgnc_symbol","wikigene_name","uniprot_genename")[ which( c("hgnc_symbol","wikigene_name","uniprot_genename") %in% listAttributes(mart)[,1] ) ][1]
      suppressWarnings({  biomart.table <-
        getBM(c(env$preferences$database.id.type, query),
              env$preferences$database.id.type,
              rownames(env$indata)[seq(1,nrow(env$indata),length.out=100)],
              mart, checkFilters=FALSE)  })
    }, silent=TRUE)

    if (is.null(biomart.table) || nrow(biomart.table) == 0)
    {
      util.warn("Invalid annotation parameters. Trying autodetection...")
      env$preferences$database.dataset <- "auto"
    }
  }

  if (env$preferences$database.dataset == "auto")
  {
    env <- pipeline.detectEnsemblDataset(env)
  }
  if (env$preferences$database.dataset == "" || env$preferences$database.id.type == "")
  {
    util.warn("Could not find valid annotation parameters.")
    util.warn("Disabling geneset & PSF analyses.")
    env$preferences$activated.modules$geneset.analysis <- FALSE
    env$preferences$activated.modules$psf.analysis <- FALSE
    return(env)
  }
  

  mart <- useMart(biomart=env$preferences$database.biomart, host=env$preferences$database.host)
  mart <- useDataset(env$preferences$database.dataset, mart=mart)

  query = c("wikigene_name","hgnc_symbol","uniprot_genename")[ which( c("wikigene_name","hgnc_symbol","uniprot_genename") %in% listAttributes(mart)[,1] ) ][1]
  
  biomart.table <- NULL
  try({
    biomart.table <- getBM(c(env$preferences$database.id.type, query,
        "description","ensembl_gene_id","chromosome_name","band","start_position"),
        env$preferences$database.id.type, rownames(env$indata), mart, checkFilters=FALSE)  

    biomart.table <- biomart.table[ which(biomart.table[,1]%in%rownames(env$indata)), ]
  }, silent=TRUE)  
  
  if (is.null(biomart.table) || nrow(biomart.table) == 0)
  {
    util.warn("Could not resolve rownames. Possibly wrong database.id.type")
    util.warn("Disabling geneset & PSF analyses.")
    env$preferences$activated.modules$geneset.analysis <- FALSE
    env$preferences$activated.modules$psf.analysis <- FALSE
  } else
  {
    h <- biomart.table[,query]
    names(h) <- biomart.table[,env$preferences$database.id.type]
    h <- h[ which(h!="") ]
    h <- h[ which(!duplicated(names(h))) ]
    env$gene.info$names[names(h)] <- h

    h <- biomart.table[,"description"]
    names(h) <- biomart.table[,env$preferences$database.id.type]
    h <- h[ which(h!="") ]
    h <- h[ which(!duplicated(names(h))) ]
    h <- sub("\\[.*\\]","",h)
    env$gene.info$descriptions[names(h)] <- h

    h <- biomart.table[,"ensembl_gene_id"]
    names(h) <- biomart.table[,env$preferences$database.id.type]
    h <- h[ which(h!="") ]
    h <- h[ which(!duplicated(names(h))) ]
    env$gene.info$ids[names(h)] <- h

    h <- biomart.table[,"chromosome_name"]
    names(h) <- biomart.table[,env$preferences$database.id.type]
    h <- h[ which(h!="") ]
    h <- h[ which(!duplicated(names(h))) ]
    env$gene.info$chr.name[names(h)] <- h
    
    h <- gsub("\\..*$","", biomart.table[,"band"])
    names(h) <- biomart.table[,env$preferences$database.id.type]
    h <- h[ which(h!="") ]
    h <- h[ which(!duplicated(names(h))) ]
    env$gene.info$chr.band[names(h)] <- h
    
    h <- biomart.table[,"start_position"]
    names(h) <- biomart.table[,env$preferences$database.id.type]
    h <- h[ which(h!="") ]
    h <- h[ which(!duplicated(names(h))) ]
    env$gene.info$chr.start[names(h)] <- h
    
    env$gene.info$ensembl.mapping <- unique(biomart.table[,c(env$preferences$database.id.type,"ensembl_gene_id")])

    gene.positions.table <- cbind( env$gene.info$chr.name, env$gene.info$chr.band )
    gene.positions.table <- gene.positions.table[ which( gene.positions.table[,1] != "" & gene.positions.table[,2] != "" ) , ]
    skip.chrnames <- names(which(table(env$gene.info$chr.name) < 20)) # filter low abundand chromosome information
    skip.chrnames <- union( skip.chrnames,  unique(env$gene.info$chr.name)[ which(nchar( unique(env$gene.info$chr.name) ) >= 8) ]   )  # filter strange CHR_ chromosomes
    gene.positions.table <- gene.positions.table[ which( !gene.positions.table[,1] %in% skip.chrnames ) , ]
      
    if(nrow(gene.positions.table)>0)
    {
      env$chromosome.list <- tapply(rownames(gene.positions.table), gene.positions.table[,1], c)
      env$chromosome.list <- lapply(env$chromosome.list, function(x) { tapply(x, gene.positions.table[x,2], c) })
    }

  }
  
  if( env$preferences$database.dataset!="hsapiens_gene_ensembl")
  {
    util.warn("Disabling PSF analysis (available for human data only).")
    env$preferences$activated.modules$psf.analysis <- FALSE 
  }

  if (!env$preferences$activated.modules$geneset.analysis) {
    return(env)
  }
  
  suppressWarnings({  biomart.table <- getBM(c("ensembl_gene_id", "go_id", "name_1006", "namespace_1003"), "ensembl_gene_id", unique(env$gene.info$ensembl.mapping$ensembl_gene_id), mart, checkFilters=FALSE) })
  biomart.table <- biomart.table[which( apply(biomart.table,1,function(x) sum(x=="") ) == 0 ),]  
  env$gs.def.list <- tapply(biomart.table[,1], biomart.table[,2], c)
  env$gs.def.list <- lapply(env$gs.def.list, function(x) { list(Genes=x, Type="") })


  if (length(env$gs.def.list) > 0)
  {
    ## simple small-gs-filtering
    env$gs.def.list <- env$gs.def.list[ which( sapply(env$gs.def.list, function(x) length(x$Genes)) >= 20 ) ]

    biomart.table[,4] <- sub("biological_process", "BP", biomart.table[,4])
    biomart.table[,4] <- sub("molecular_function", "MF", biomart.table[,4])
    biomart.table[,4] <- sub("cellular_component", "CC", biomart.table[,4])

    for (i in seq_along(env$gs.def.list))
    {
      o <- match(names(env$gs.def.list)[i], biomart.table[,2])
      names(env$gs.def.list)[i] <- biomart.table[o, 3]
      env$gs.def.list[[i]]$Type <- biomart.table[o, 4]
    }

    env$gs.def.list <- env$gs.def.list[order(names(env$gs.def.list))]
  } else
  {
    util.warn("No GO annotation found")
  }

  if(length(env$chromosome.list)>0)
  {
    chr.gs.list <- lapply(env$chromosome.list, function(x) { list(Genes=env$gene.info$ids[unlist(x)], Type="Chr") })
    names(chr.gs.list) <- paste("Chr", names(env$chromosome.list))
    env$gs.def.list <- c(env$gs.def.list, chr.gs.list)
  }


  # load custom genesets
  data(opossom.genesets)
  env$gs.def.list <- c(env$gs.def.list, opossom.genesets)

  env$gs.def.list <- lapply(env$gs.def.list, function(x) {
    x$Genes <- intersect(x$Genes, unique(env$gene.info$ensembl.mapping$ensembl_gene_id))
    return(x)
  })


  env$gs.def.list <- env$gs.def.list[ which(sapply(sapply(env$gs.def.list, head, 1), length) >= 10) ]
  env$gs.def.list <- env$gs.def.list[ which(sapply(names(env$gs.def.list),nchar) < 60 ) ]

  if (length(env$gs.def.list) > 0)
  {
    env$gs.def.list <- lapply(env$gs.def.list, function(x) { names(x$Genes) = NULL; return(x) })
    util.info(length(env$gs.def.list), "annotated gene sets derived for functional analysis")
  } else
  {
    env$preferences$activated.modules$geneset.analysis <- FALSE
    util.warn("No Geneset information -> turning off GS analysis")
  }
  return(env)
}
