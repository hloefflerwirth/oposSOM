
indata = samples.GSZ.scores
rm( indata.original )
rm( colramp )


preferences$dataset.name = paste( preferences$dataset.name, "GenesetSOM" )

preferences$error.model = "all.samples"

preferences$ensembl.dataset = ""
preferences$ensembl.rowname.ids = ""
preferences$geneset.analysis = F

preferences$feature.mean.normalization = F
preferences$training.extension = 4


source("R/source/pipeline.r")









#require.bioconductor("GO.db")

#unlist( as.list(GOBPOFFSPRING["GO:0033059"]) )

#gs.def.list = list()

#as.list(GOTERM["GO:0033059"])





