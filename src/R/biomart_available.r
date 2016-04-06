biomart.available <- function()
{
  test.table <- NULL

  test.rownames <- c("ENSG00000115415",
                     "ENSG00000049541",
                     "ENSG00000143384",
                     "ENSG00000215301",
                     "ENSG00000111348",
                     "ENSG00000114978",
                     "ENSG00000166888",
                     "ENSG00000154473",
                     "ENSG00000123136",
                     "ENSG00000162511")

  try({
    mart <- useMart(biomart=preferences$database.biomart, host=preferences$database.host)
    mart <- useDataset("hsapiens_gene_ensembl", mart=mart)

    suppressWarnings({  test.table <- getBM(c("hgnc_symbol", "ensembl_gene_id"),
                        "ensembl_gene_id",
                        test.rownames,
                        mart=mart)  })
  }, silent=TRUE)

  return(is.data.frame(test.table) && nrow(test.table) > 0)
}
