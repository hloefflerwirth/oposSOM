

gene.names = rep( "", nrow( indata ) )
names( gene.names ) = rownames(indata)

gene.descriptions = rep( "", nrow( indata ) )
names( gene.descriptions ) = rownames(indata)

gene.ids = rep( "", nrow( indata ) )
names( gene.ids ) = rownames(indata)

gene.positions = rep( "", nrow( indata ) )
names( gene.positions ) = rownames(indata)
gene.positions.table = matrix(0,0,0)
gene.positions.list = list()


if( preferences$ensembl.dataset != "" && preferences$ensembl.rowname.ids != "" )
{

	require.bioconductor( "biomaRt" )



	mart<-useMart('ensembl')
	mart<-useDataset(preferences$ensembl.dataset,mart=mart)

	biomart.table = getBM( c( preferences$ensembl.rowname.ids, "external_gene_id", "description", "ensembl_gene_id", "chromosome_name","band" ) , preferences$ensembl.rowname.ids, rownames(indata), mart, checkFilters=F )
	

	if( nrow(biomart.table) == 0 )
	{
		preferences$geneset.analysis = F
		cat( "!!Could not resolve rownames. Possibly wrong ensembl.rowname.ids!!\n" ); flush.console()

	} else
	{		
		h = biomart.table[,2]
		names(h) = biomart.table[,1]
		gene.names[ as.character( unique(biomart.table[,1]) ) ] = h[ as.character( unique(biomart.table[,1]) ) ]
		
		h = biomart.table[,3]
		names(h) = biomart.table[,1]
		gene.descriptions[ as.character( unique(biomart.table[,1]) ) ] = h[ as.character( unique(biomart.table[,1]) ) ]
	
		h = biomart.table[,4]
		names(h) = biomart.table[,1]
		gene.ids[ as.character( unique(biomart.table[,1]) ) ] = h[ as.character( unique(biomart.table[,1]) ) ]
		gene.ids = gene.ids[ which( gene.ids != "" ) ]
		gene.ids = gene.ids[ which( names(gene.ids) %in% rownames(indata) ) ] 
	
		h = paste( biomart.table[,"chromosome_name"], gsub("\\..*$","", biomart.table[,"band"] ) )
		names(h) = biomart.table[,1]
		gene.positions[ as.character( unique(biomart.table[,1]) ) ] = h[ as.character( unique(biomart.table[,1]) ) ]
		gene.positions = gene.positions[ which( gene.positions != "" ) ]
	
		gene.positions.table = do.call( rbind, strsplit( gene.positions, " " ) )

			
		junk.chrnames = names( which(  table(gene.positions.table[,1]) < 20 ) )
		junk.chrnames = union( junk.chrnames, names( which( tapply( gene.positions.table[,2], gene.positions.table[,1], function(x) length(unique(x)) ) == 1 ) ) )
		gene.positions = gene.positions[ which( ! gene.positions.table[,1] %in% junk.chrnames ) ]		
		gene.positions.table = gene.positions.table[ which( ! gene.positions.table[,1] %in% junk.chrnames ), ]	
		
		
		gene.positions.list = tapply( rownames(gene.positions.table), gene.positions.table[,1], c )
		gene.positions.list = lapply( gene.positions.list, function(x){ tapply( x, gene.positions.table[x,2], c )  })
	
	}






	if( preferences$geneset.analysis )
	{
#		require.bioconductor( "GO.db" )
		

		unique.protein.ids = unique(gene.ids)



		biomart.table = getBM( c( "ensembl_gene_id", "go_id", "name_1006", "namespace_1003" ) , "ensembl_gene_id", unique.protein.ids, mart, checkFilters=F )
		gs.def.list = tapply( biomart.table[,1], biomart.table[,2], c  )
		gs.def.list = lapply( gs.def.list, function(x){ list( Genes=x, Type="" ) } )


		gs.def.list = gs.def.list[ - which( names(gs.def.list) == "" ) ]


		if( length(gs.def.list) > 0 )
		{

	                                             	
			cat( "Download of", length(gs.def.list), "GO sets with", sum( sapply( sapply( gs.def.list, head, 1 ), length ) ), "entries\n" ); flush.console()
	
	
			
			## Reduce number of GO sets
	

		
# 			# Construct GO list with sets only containing regulated genes
# 	
# 			sig = na.omit( gene.ids[  names( which( apply( abs(indata), 1, max ) >= quantile( apply( abs(indata), 1, max ), 0.9 ) ) )  ] )
# 	
# 			for( i in 1:length(gs.def.list) )
# 			{
# 				gs.def.list[[i]]$Genes = gs.def.list[[i]]$Genes[  which( gs.def.list[[i]]$Genes %in% sig )  ]
# 			}
# 			gs.def.list = gs.def.list[ - which( sapply( sapply( gs.def.list, head, 1 ), length ) < 1 ) ]
# 	
# 	
# 	
# 	
# 	
# 			# Construct GO tree in igraph
# 	
# 			GO.parent.list = c( as.list(GOBPPARENTS), as.list(GOMFPARENTS), as.list(GOCCPARENTS) )
# 			GO.parent.list = GO.parent.list[  names(GO.parent.list) %in% names(gs.def.list)  ]
# 	
# 			for( i in 1:length(GO.parent.list) )
# 			{
# 				GO.parent.list[[i]] = GO.parent.list[[i]][  which( GO.parent.list[[i]] %in% names(gs.def.list) )  ]
# 			}
# 			GO.parent.list = GO.parent.list[ - which( sapply( GO.parent.list, length ) < 1 ) ]
# 	
# 	
# 	
# 			edge.matrix = matrix( 0, 0, 2 )
# 			for( i in 1:length(GO.parent.list) )
# 			{	
# 				for( j in 1:length(GO.parent.list[[i]]) )
# 				{
# 					edge.matrix = rbind( edge.matrix,    c( names(GO.parent.list)[i], GO.parent.list[[i]][j] )    )
# 				}
# 			}
# 			g = graph.edgelist( edge.matrix )
# 	
# 	
# 	
# 	
# 			## Elim algorithm to reduce redundancy		
# 	
# 			level.nodes = V(g)[  which( degree( g, mode="in" ) == 0 ) - 1  ]
# 			remove.list = c()
# 			while( length(level.nodes) > 0 )
# 			{
# 				next.remove.list = c()
# 	
# 				for( i in 1:length( gs.def.list[ level.nodes ] ) )
# 				{
# 					if( length( which( gs.def.list[ level.nodes ][[i]]$Genes %in% remove.list) ) )
# 					{
# 						gs.def.list[ level.nodes ][[i]]$Genes = gs.def.list[ level.nodes ][[i]]$Genes[ -which( gs.def.list[ level.nodes ][[i]]$Genes %in% remove.list) ]
# 					}
# 	
# 					next.remove.list = c( next.remove.list, gs.def.list[ level.nodes ][[i]]$Genes )
# 				}
# 				remove.list = next.remove.list	
# 				level.nodes = setdiff(  unique( unlist( neighborhood( g, order=1, nodes=level.nodes, mode="out" ) ) ),   level.nodes )
# 			}
# 			gs.def.list = gs.def.list[ - which( sapply( sapply( gs.def.list, head, 1 ), length ) < 1 ) ]
# 	
# 	
# 
# 
# 			## Used reduced GO tree to filter sets for analysis	
# 
# 			name.list = names(gs.def.list)
# 
# 			gs.def.list = tapply( biomart.table[,1], biomart.table[,2], c  )
# 			gs.def.list = lapply( gs.def.list, function(x){ list( Genes=x, Type="" ) } )
# 			gs.def.list = gs.def.list[ name.list ]
# 
# 			gs.def.list = gs.def.list[ - which( sapply( sapply( gs.def.list, head, 1 ), length ) < 5 ) ]
	
			
	
	
			## simple small-gs-filtering
		
			gs.def.list = gs.def.list[ - which( sapply( sapply( gs.def.list, head, 1 ), length ) < 10 ) ]

	
	
	
	
			cat( "Filtered to", length(gs.def.list), "GO sets with", sum( sapply( sapply( gs.def.list, head, 1 ), length ) ), "entries\n" ); flush.console()
	

	
	
			biomart.table[,4] = sub( "biological_process","BP", biomart.table[,4] )
			biomart.table[,4] = sub( "molecular_function","MF", biomart.table[,4] )
			biomart.table[,4] = sub( "cellular_component","CC", biomart.table[,4] )
	
			for( i in 1:length(gs.def.list) )
			{
				o = match( names(gs.def.list)[i], biomart.table[,2] )		
	
				names(gs.def.list)[i] = biomart.table[ o, 3 ]
				gs.def.list[[i]]$Type = biomart.table[ o, 4 ]
			}
	
	
	
	
			gs.def.list = gs.def.list[ order( names( gs.def.list ) ) ]
	

		} else
		{
			cat( "!!No GO annotation found!!\n" ); flush.console()
		}

		
		
		
		
# 		for( i in 1:length(gene.positions.list) )
# 		{
# 			chr.gs.list = list( list( Genes=gene.ids[ unlist( gene.positions.list[[i]] ) ], Type="Chr") )
# 			chr.gs.list = c( chr.gs.list, lapply( gene.positions.list[[i]], function(x) list( Genes=gene.ids[x], Type="Chr" ) ) )
# 			names(chr.gs.list) = paste( "Chr",names(gene.positions.list)[i], names(chr.gs.list) )
# 			
# 			gs.def.list = c( gs.def.list, chr.gs.list )
# 		}
		
		chr.gs.list = lapply( gene.positions.list, function(x) list( Genes=gene.ids[unlist(x)], Type="Chr" ) )
		names(chr.gs.list) = paste( "Chr", names(gene.positions.list) )		
		gs.def.list = c( gs.def.list, chr.gs.list )

		
		
		small.gs = which( sapply( sapply( gs.def.list, head, 1 ), length ) < 10 )
		if( length(small.gs) > 0 )
		{
			gs.def.list = gs.def.list[ - small.gs ]
		}

		
		
		if( length( preferences$geneset.custom.list ) > 0 ) 
		{				
			gs.def.list = c( gs.def.list, preferences$geneset.custom.list )
		}



		gs.def.list = lapply( gs.def.list, function(x)
		{ 
			x$Genes = intersect( x$Genes, unique.protein.ids )	
			return(x)
		} )

		
		small.gs = which( sapply( sapply( gs.def.list, head, 1 ), length ) < 2 )
		if( length(small.gs) > 0 )
		{
			gs.def.list = gs.def.list[ - small.gs ]
		}



		if( length(gs.def.list) > 0 )
		{           
			gs.def.list = lapply( gs.def.list, function(x){ names(x$Genes) = NULL; return(x) } )
			gs.def.list.categories = sapply(gs.def.list,function(x){ x$Type } )
			cat( "In total",length(gs.def.list),"gene sets to be considered in analysis\n" ); flush.console()
		} else
		{
			preferences$geneset.analysis = F
			cat( "!!No Geneset information -> switch off GS analysis!!\n" ); flush.console()
		}		
		
		
	} 




} else
{

	preferences$geneset.analysis = F

}



