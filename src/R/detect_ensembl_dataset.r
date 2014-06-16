pipeline.detectEnsemblDataset <- function()
{
  preferences$database.dataset <<- ""
  util.info("Autodetecting annotation parameters")

  if (!suppressWarnings(any(is.na(as.numeric(rownames(indata))))) ||
      !suppressWarnings(any(is.na(as.numeric(colnames(indata))))))
  {
    return()
  }

  auto.datasets <-
    list("hsapiens_gene_ensembl"=c("affy_hg_u133_plus_2",
                                   "affy_hg_u133a",
                                   "affy_hugene_1_0_st_v1",
                                   "ipi",
                                   "refseq_mrna"),
         "mmusculus_gene_ensembl"=c("affy_moe430a",
                                    "affy_mogene_1_0_st_v1",
                                    "illumina_mousewg_6_v2"),
         "rnorvegicus_gene_ensembl"=c("affy_rae230a"),
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
        biomart.table <-
          getBM(c(id, "external_gene_id"), id,
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
