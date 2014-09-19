

csv.file.to.compare = "oposSOM Geneset Collection.csv"



library(oposSOM)


data(opossom.genesets)
opossom.genesets = names(opossom.genesets)
opossom.genesets = opossom.genesets[ which( opossom.genesets != "" ) ]




new.genesets = read.csv2(csv.file.to.compare, as.is=T)
new.genesets = new.genesets$Name
new.genesets = new.genesets[ which( new.genesets != "" ) ]




cat("Genesets in oposSOM, but not in your collection:\n" )
cat( paste( setdiff( opossom.genesets, new.genesets ), collapse="\n" ), "\n" )

cat("Genesets in your collection, but not in oposSOM:\n" )
cat( paste( setdiff( new.genesets, opossom.genesets ), collapse="\n" ), "\n" )






