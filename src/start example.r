
  rm(list=ls(all=TRUE))

  # Biomart parameter: #
  ######################

  # Human:
    # dataset: "hsapiens_gene_ensembl"
      # "affy_hg_u133a"
      # "affy_hg_u133_plus_2"
      # "affy_hugene_1_0_st_v1"
      # "ipi"
      # "refseq_mrna" eg "NM_173201"

  # Mouse:
    # dataset: "mmusculus_gene_ensembl"
    # rowname.ids:
      # "affy_moe430a"
      # "affy_mogene_1_0_st_v1"
      # "illumina_mousewg_6_v2"

  # Rat:
    # dataset: "rnorvegicus_gene_ensembl"
    # rowname.ids:
      # "affy_rae230a"

  # Yeast:
    # dataset: "scerevisiae_gene_ensembl"
    # rowname.ids:
      # "affy_yg_s98"
      # "affy_yeast_2"


  # common rowname.ids:
      # "ensembl_gene_id'
      # "entrezgene"
      # "hgnc_symbol"


  # Parameters

  pipeline <- opossom.new(list(dataset.name        = "Test",
                               error.model         = "all.samples",  # "replicates", "all.samples" or "groups"
                               dim.som1            = 20,
                               dim.som2            = 20,
                               training.extension  = 1,

                               rotate.som1         = 0,
                               flip.som1           = F,

                               ensembl.dataset      = "hsapiens_gene_ensembl",
                               ensembl.rowname.ids = "affy_hg_u133_plus_2",

                               geneset.analysis     = T,
                               geneset.analysis.exact = F,

                               max.parallel.cores = 4,
                               sample.spot.cutoff  = 0.65,
                               summary.spot.core  = 5,
                               summary.spot.threshold = 0.95,
                               group.spot.core = 5,
                               group.spot.threshold = 0.75,

                               feature.mean.normalization = T,
                               sample.quantile.normalization = T ))

#                differences.list = list("group1 vs group2" = list( c(1:3), c(4:6) ),
#                                         "group1 vs group3" = list( c(1:3), 7 ) ) )






  # Load Data

  pipeline$indata <- .......

    # new hook readin:

    # indata = log10( read.table( "allExpresssionMeasures.dat", sep="\t", comment="", header=T, row.names=1, check.names=F) )
    # colnames(indata) = sub( "/blabla/", "", colnames(indata) )
    # colnames(indata) = sub( ".CEL", "", colnames(indata) )



  # Set group labels and colors for each column in data

  pipeline$group.labels <- .....

    # group.colors = c("col1","col2",...)[ match( group.labels, unique(group.labels) ) ]



  # optional for hook-preprocessed data: quantile + camel analysis

    #    source("R/quantile_normalization.r")
    #    indata = Quantile.Normalization( indata )

    #   source("R/camel_analysis.r")
    #   train.camel = Camel.Analysis( indata )
    #   indata = train.camel$Corrected.Data * train.camel$Percent.S
    #   indata.original = train.camel$Corrected.Data
    #   rm( train.camel )



  # Lets rock!

  opossom.run(pipeline)




