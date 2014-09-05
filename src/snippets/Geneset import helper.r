


library(biomaRt)

mart<-useMart('ensembl')
mart<-useDataset('hsapiens_gene_ensembl',mart=mart)



grep("NM_", listFilters(mart)[,2], value=T )






tab = read.table( "lala.csv" )
ids = as.character( tab[,1] )
ids = ids[ which(ids!="") ]



ids = c("","","","","","","","","","","","","" )

ids = c("","","","","","","","","","","","","" )







biomart.table = getBM( c("ensembl_gene_id", "hgnc_symbol"), "hgnc_symbol", ids, mart, checkFilters=F )
biomart.table = getBM( c("ensembl_gene_id", "affy_hg_u133_plus_2"), "affy_hg_u133_plus_2", ids, mart, checkFilters=F )
biomart.table = getBM( c("ensembl_gene_id", "affy_hugenefl"), "affy_hugenefl", ids, mart, checkFilters=F )
biomart.table = getBM( c("ensembl_gene_id", "refseq_mrna"), "refseq_mrna", ids, mart, checkFilters=F )


biomart.table
ids[ which( ! ids %in% biomart.table[,2] ) ]



write.csv2( paste( unique( biomart.table[,1] ), collapse="," ), file="GS.csv" )


















sets = list(

	
  CARO_OxPhos_in_DLBCL_UP = c( "ACADVL","ACAT1","AK2","ALDH18A1","ALDH1B1","ALDH2","ALDH5A1",
                        "ATP5A1","ATP5B","ATP5C1","ATP5H","ATP5I","ATP5J","ATP5L","ATP5O",
                        "CPT2","ECHS1","ETFA","ETFB","HADHA","HADHB","IDH2","IDH3A",
                        "NDUFA5","NDUFA8","NDUFA9","NDUFAF2","NDUFS1","NDUFS2","NDUFS3",
                        "NDUFV1","NDUFV2","NDUFV3","OAT","OGDH","OXCT1","PDHB","SDHA",
                        "SDHB","SHMT2","SOD2","VDAC3" ),
  
  CARO_OxPhos_vs_BCR_UP = c( "PDHB","NDUFA5","NDUFS3","ACAT1","ETFA","MTHFD2","SOD2" )  

)



sets = lapply( sets, function(x)
{

	biomart.table = getBM( c("ensembl_gene_id", "hgnc_symbol"), "hgnc_symbol", x, mart, checkFilters=F )
	
	ids = unique( biomart.table[,1] )
	ids = ids[ grep( "ENS", ids ) ]
	
	return(  paste( ids, collapse="," )  )
	
})





write.csv2( do.call(rbind,sets), file="GS.csv" )




