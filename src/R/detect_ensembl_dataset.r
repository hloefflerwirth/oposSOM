pipeline.detectEnsemblDataset <- function()
{
  preferences$database.dataset <<- ""
  util.info("Autodetecting annotation parameters")

  auto.datasets <-
    list("hsapiens_gene_ensembl"=c("affy_hg_u133_plus_2",
                                   "affy_hugene_1_0_st_v1",
                                   "affy_hugene_2_0_st_v1",
                                   "ipi",
                                   "refseq_mrna"),
         "mmusculus_gene_ensembl"=c("affy_moe430a",
                                    "affy_mogene_1_0_st_v1",
                                    "affy_mogene_2_1_st_v1",
                                    "illumina_mousewg_6_v2"),
         "rnorvegicus_gene_ensembl"=c("affy_rae230a",
                                      "affy_ragene_1_0_st_v1",
                                      "affy_ragene_2_0_st_v1"),
         "scerevisiae_gene_ensembl"=c("affy_yg_s98",
                                      "affy_yeast_2"))

  auto.rowname.ids <- c("ensembl_gene_id", "entrezgene", "hgnc_symbol")

  mart <- useMart('ensembl')

  for (ds in names(auto.datasets))
  {
    mart <- useDataset(ds, mart=mart)

    for (id in c(auto.datasets[[ds]], auto.rowname.ids))
    {
      try({
        query = c("wikigene_name","hgnc_symbol","uniprot_genename")[ which( c("wikigene_name","hgnc_symbol","uniprot_genename") %in% listAttributes(mart)[,1] ) ][1]
        biomart.table <-
          getBM(c(id, query), id,
                rownames(indata)[seq(1,nrow(indata),length.out=100)],
                mart, checkFilters=FALSE)

        if (nrow(biomart.table) > 0)
        {
          util.info("Detected annotation dataset:", ds)
          util.info("Detected annotation filter:", id)
          preferences$database.dataset <<- ds
          preferences$database.id.type <<- id
          return()
        }
      }, silent=TRUE)
    }
  }
}
