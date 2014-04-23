

  geneset.custom.matrix = read.csv2("data/geneset_collection.csv", as.is=T )


  # collapse gene sets in more than one cell
  for( i in 1:nrow(geneset.custom.matrix) )
  {

    n.id.cells = sum( sapply( geneset.custom.matrix[i,c(4:ncol(geneset.custom.matrix))], nchar ) > 0 )
    if( n.id.cells > 1)
    {
      geneset.custom.matrix[i,4] = paste( geneset.custom.matrix[i,c(4:(4+n.id.cells-1))], collapse="," )
    }

  }



  selection = which( geneset.custom.matrix[,"Type"] %in% geneset.custom.selection )



  geneset.custom.list = apply( geneset.custom.matrix[ selection, ], 1, function(x){ list( Genes=as.vector(x["Genes"]), Type=as.vector(x["Type"]) ) })

  names(geneset.custom.list) = geneset.custom.matrix[ selection, "Name" ]



  geneset.custom.list = lapply( geneset.custom.list, function(x)
  {
    x$Genes = strsplit( x$Genes, "," )[[1]]
    return(x)
  } )


  rm(geneset.custom.matrix)
