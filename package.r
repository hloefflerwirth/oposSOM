
fast <- TRUE
install <- TRUE

# Package directories
dirname <- "oposSOM"
datadir <- file.path(dirname, "data")

# Remove existing package dir
if (file.exists(dirname)) {
  cat("* remove existing package dir:", dirname, "\n")
  unlink(dirname, recursive=T)
}
dir.create(dirname)

# Copy everything we need
for (f in c("vignettes", "man", "inst", "DESCRIPTION", "NAMESPACE","NEWS")) {
  cat("* copy", f, "\n")
  file.copy(file.path("src", f), file.path(dirname), recursive=T)
}

cat("* copy R\n")

dir.create(paste(dirname,"/R",sep=""))
dummy <- file.copy( paste( "src/R/",dir( "src/R" ),sep=""), paste(dirname,"/R/",sep="") )


# Move data into the package
cat("* create data dir:", datadir, "\n")
dir.create(datadir)

filename <- "opossom.tissues.RData"
cat("* copy", filename, "\n")
file.copy(file.path("src", "data", filename), datadir) -> x


### build gene set collection
  
  ### from csv collection
  filename <- file.path("src", "data", "geneset_collection.csv")
  cat("* load", filename, "\n")
  genesets <- read.csv2(filename, as.is=T)
  
  # collapse gene sets in more than one cell
  for(i in 1:nrow(genesets)) {
    cells <- sum(sapply(genesets[i, c(4:ncol(genesets))], nchar) > 0)
  
    if (cells > 1) {
      genesets[i,4] = paste(genesets[i, c(4:(4+cells-1))], collapse=",")
    }
  }
  
  opossom.genesets <- apply(genesets, 1, function(x) {
    list(Genes=as.vector(x["Genes"]), Type=as.vector(x["Type"]))
  })
  names(opossom.genesets) <- genesets[,"Name"]
  
  opossom.genesets <- lapply(opossom.genesets, function(x) {
    x$Genes = strsplit(x$Genes, ",")[[1]]
    return(x)
  })

  opossom.genesets <- opossom.genesets[ which( names(opossom.genesets) != "" ) ]

  

  ### from csv collection - inhouse
  filename <- file.path("src", "data", "geneset_collection - inhouse.csv")
  cat("* load", filename, "\n")
  genesets <- read.csv2(filename, as.is=T)
  
  # collapse gene sets in more than one cell
  for(i in 1:nrow(genesets)) {
    cells <- sum(sapply(genesets[i, c(4:ncol(genesets))], nchar) > 0)
    
    if (cells > 1) {
      genesets[i,4] = paste(genesets[i, c(4:(4+cells-1))], collapse=",")
    }
  }
  
  n = genesets$Name
  genesets <- apply(genesets, 1, function(x) {
    list(Genes=as.vector(x["Genes"]), Type=as.vector(x["Type"]))
  })
  names(genesets) <- n
  
  genesets <- lapply(genesets, function(x) {
    x$Genes = strsplit(x$Genes, ",")[[1]]
    return(x)
  })
  
  genesets <- genesets[ which( names(genesets) != "" ) ]  
  
  
  opossom.genesets <- c( opossom.genesets, genesets )
  
  
  
  

  ### from gsea files
  library(biomaRt)
  mart<-useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")


  filename <- file.path("src", "data", "h.all.v5.0.entrez.gmt")
  cat("* load", filename, "\n")
  genesets <- readLines(filename)
  genesets <- sapply( lapply( genesets, strsplit, "\t" ), head, 1 )
  names(genesets) <- sapply(genesets, head,1)
  genesets  <- sapply(genesets, function(x) x[-c(1,2)] )
  
  biomart.table <- getBM( c( "entrezgene", "ensembl_gene_id" ) , "entrezgene", unique(unlist(genesets)), mart )
  genesets <- sapply(genesets, function(x) biomart.table[match(x,biomart.table[,1]),2] )
  genesets <- lapply(genesets, function(x) as.vector(na.omit(x)) )
  genesets <- genesets[ which( sapply(genesets,length) > 0 ) ]
  genesets <- lapply(genesets, function(x) list(Genes=x, Type="HM") )

  opossom.genesets <- c( opossom.genesets, genesets )
  
  
  filename <- file.path("src", "data", "c2.all.v5.0.entrez.gmt")
  cat("* load", filename, "\n")
  genesets <- readLines(filename)
  genesets <- sapply( lapply( genesets, strsplit, "\t" ), head, 1 )
  names(genesets) <- sapply(genesets, head,1)
  genesets  <- sapply(genesets, function(x) x[-c(1,2)] )
  
  biomart.table <- getBM( c( "entrezgene", "ensembl_gene_id" ) , "entrezgene", unique(unlist(genesets)), mart )
  genesets <- sapply(genesets, function(x) biomart.table[match(x,biomart.table[,1]),2] )
  genesets <- lapply(genesets, function(x) as.vector(na.omit(x)) )
  genesets <- genesets[ which( sapply(genesets,length) > 0 ) ]
  genesets <- lapply(genesets, function(x) list(Genes=x, Type="GSEA C2") )

  opossom.genesets <- c( opossom.genesets, genesets )
  


filename <- file.path(datadir, "opossom.genesets.RData")
cat("* save", filename, "\n")
save(opossom.genesets, file=filename, compress="xz")



# Build and Check
cat("\n>> R CMD build --resave-data", dirname, "\n")
system(paste("R CMD build --resave-data", dirname))

unlink("oposSOM",recursive=T)



if (!fast) {
  cat(">> R CMD check", dirname, "\n")
  system(paste("R CMD check", dirname))
}

# Install
if (install) {
	# remove.packages("oposSOM")
  install.packages( paste( "H:\\Eigene Dateien\\Entwicklung\\R Scripts\\SOM Pipeline\\oposSOM\\Package Build\\",
												  	sort( grep( "oposSOM_", dir(), value=TRUE), decreasing=TRUE )[1], sep="" ),
														repos=NULL, type="source" )
}



system("timeout /T 600")



