

# 	bootstrap.col = rgb2hsv( col2rgb(group.colors) )
# 	bootstrap.col[2,] = bootstrap.col[2,] - ( 0.8 * as.numeric( group.bootstrap.score < 80 ) + 0.19 * as.numeric( group.bootstrap.score < 50 ) )
# 	bootstrap.col[ which(bootstrap.col < 0 ) ] = 0
# 	bootstrap.col = hsv( bootstrap.col[1,], bootstrap.col[2,], bootstrap.col[3,] )

	bootstrap.col = group.colors
	bootstrap.col[ which(group.bootstrap.score < 80) ] = apply( ( col2rgb(bootstrap.col[ which(group.bootstrap.score < 80) ]) + 0.7 * ( 255 - col2rgb(bootstrap.col[ which(group.bootstrap.score < 80) ]) ) ) /255, 2, function(x) rgb(x[1],x[2],x[3]) )
	bootstrap.col[ which(group.bootstrap.score < 50) ] = apply( ( col2rgb(bootstrap.col[ which(group.bootstrap.score < 50) ]) + 0.9 * ( 255 - col2rgb(bootstrap.col[ which(group.bootstrap.score < 50) ]) ) ) /255, 2, function(x) rgb(x[1],x[2],x[3]) )

	
	
	
	pdf( paste( files.name, "- Results/2nd lvl Metagene Analysis/Correlation Analysis.pdf" ), 29.7/2.54, 21/2.54 )


	for( i in 1:length(metagene.filter.list) )
	{	
		
		s = metadata[ metagene.filter.list[[i]]$s , ]
		par( mar=c(1,1,1,1) )


		
		
		
		# Maximum Correlation Tree
		

		adj.matrix = cor( s ) * -1
		g = graph.adjacency( adj.matrix, weighted=T,  mode="undirected" )
		stg = minimum.spanning.tree( g ) 

		layout=layout.kamada.kawai(stg)

		plot( stg, layout=layout, vertex.size=5, vertex.label = colnames(indata), vertex.label.cex=if(ncol(indata)<100) 1.2 else 0.6, vertex.color=group.colors, main=metagene.filter.list[[i]]$n )
			legend( "bottomright", as.character(unique( group.labels )), cex=0.5, text.col=unique.group.colors, bg="white" )	
			box()

		plot( stg, layout=layout, vertex.size=5, vertex.label = rep("",ncol(indata)) , vertex.label.cex=if(ncol(indata)<100) 1.2 else 0.6, vertex.color=group.colors, main=metagene.filter.list[[i]]$n )
			legend( "bottomright", as.character(unique( group.labels )), cex=0.5, text.col=unique.group.colors, bg="white" )	
			box()

		plot( stg, layout=layout, vertex.size=5, vertex.label = rep("",ncol(indata)) , vertex.label.cex=if(ncol(indata)<100) 1.2 else 0.6, vertex.color=bootstrap.col, main=metagene.filter.list[[i]]$n )
			legend( "bottomright", as.character(unique( group.labels )), cex=0.5, text.col=unique.group.colors, bg="white" )	
			box()
		
			
		
		# collapse MCT nodes
		
# 		if( length(unique(group.labels)) > 1 )
# 		{
# 			suppressMessages({	adj.matrix = as.matrix( get.adjacency(stg) ) })			
# 	
# 			hypernode.strength = rep( 1, ncol(indata) )
# 			names(hypernode.strength) = colnames(indata)
# 	
# 			sample.i = 1
# 			while( sample.i <= ncol(indata) )
# 			{
# 				sample = colnames(indata)[sample.i]
# 				
# 				if( sample %in% colnames(adj.matrix) )
# 				{
# 					connected.samples = names( which( adj.matrix[sample,] > 0 ) )
# 				
# 					if( sum( group.labels[connected.samples] == group.labels[sample]) > 0 )
# 					{					
# 						connected.samples = connected.samples[ which( group.labels[connected.samples] == group.labels[sample] ) ]
# 						
# 						adj.matrix[sample,] = adj.matrix[sample,] + apply( adj.matrix[ connected.samples, ,drop=F], 2, sum )
# 						adj.matrix[,sample] = adj.matrix[sample,]
# 						
# 						adj.matrix = adj.matrix[ -which( rownames(adj.matrix) %in% connected.samples ), -which( colnames(adj.matrix) %in% connected.samples ) ]
# 						diag(adj.matrix) = 0
# 						
# 						hypernode.strength[sample] = hypernode.strength[sample] + length(connected.samples)
# 						hypernode.strength = hypernode.strength[ -which( names(hypernode.strength) %in% connected.samples ) ]
# 					}
# 					else
# 					{
# 						sample.i = sample.i + 1
# 					}
# 				}	
# 				else
# 				{
# 					sample.i = sample.i + 1
# 				}
# 			}
# 			
# 			
# 			
# 			g = graph.adjacency( adj.matrix, mode="undirected" )
# 			plot( g, layout=layout[match( colnames(adj.matrix), colnames(indata) ),], vertex.size=5-(sqrt(1/max(table(group.labels)))*30)+sqrt(hypernode.strength/max(table(group.labels)))*30, vertex.label = hypernode.strength , vertex.label.cex=1.2, vertex.color=group.colors[rownames(adj.matrix)], vertex.label.color="white", main=metagene.filter.list[[i]]$n )
# 				legend( "bottomright", as.character(unique( group.labels )), cex=0.5, text.col=unique.group.colors, bg="white" )	
# 				box()
# 			
# 			plot( g, layout=layout.kamada.kawai, vertex.size=5-(sqrt(1/max(table(group.labels)))*30)+sqrt(hypernode.strength/max(table(group.labels)))*30, vertex.label = hypernode.strength , vertex.label.cex=1.2, vertex.color=group.colors[rownames(adj.matrix)], vertex.label.color="white", main=metagene.filter.list[[i]]$n )
# 				legend( "bottomright", as.character(unique( group.labels )), cex=0.5, text.col=unique.group.colors, bg="white" )	
# 				box()
# 			
# 		}	
			
			

		
		# Collelation Backbone
		
		
		adj.matrix = cor( s )
		diag(adj.matrix) = 0
		
		adj.matrix = apply( adj.matrix, 1, function(x) 
		{ 		
			x[order(x,decreasing=T)[-c(1:2)]] = 0 
			return(x)
		} )
		
		adj.matrix[ which( adj.matrix < 0.5 ) ] = 0
		
		
		g = graph.adjacency( adj.matrix, weighted=T,  mode="undirected" )	
		layout=layout.fruchterman.reingold(g)
		
		plot( g, layout=layout, vertex.size=5, vertex.label = colnames(indata), vertex.label.cex=if(ncol(indata)<100) 1.2 else 0.6, vertex.color=group.colors, main=metagene.filter.list[[i]]$n )
		legend( "bottomright", as.character(unique( group.labels )), cex=0.5, text.col=unique.group.colors, bg="white" )	
		box()
		
		plot( g, layout=layout, vertex.size=5, vertex.label = rep("",ncol(indata)) , vertex.label.cex=if(ncol(indata)<100) 1.2 else 0.6, vertex.color=group.colors, main=metagene.filter.list[[i]]$n )
		legend( "bottomright", as.character(unique( group.labels )), cex=0.5, text.col=unique.group.colors, bg="white" )	
		box()
		
		plot( g, layout=layout, vertex.size=5, vertex.label = rep("",ncol(indata)) , vertex.label.cex=if(ncol(indata)<100) 1.2 else 0.6, vertex.color=bootstrap.col, main=metagene.filter.list[[i]]$n )
		legend( "bottomright", as.character(unique( group.labels )), cex=0.5, text.col=unique.group.colors, bg="white" )	
		box()		
		
		
		
		
		
		
		
		# Correlation Network

		adj.matrix = cor( s )
		diag(adj.matrix) = 0
		adj.matrix[ which( adj.matrix < 0.5 ) ] = 0

		if( max(adj.matrix) > 0 )
		{
			g = graph.adjacency( adj.matrix, weighted=T,  mode="undirected" )

#			layout=layout.fruchterman.reingold.grid( g )#, start=matrix( 1:(2*ncol(indata)),ncol=2), niter=1000 )
			layout=layout.fruchterman.reingold( g, start=matrix( 1:(2*ncol(indata)),ncol=2), niter=1000 )
			
			plot( g, layout=layout, vertex.size=ifelse( ncol(indata)<250, 5, 3 ), vertex.label = colnames(indata), vertex.label.cex=if(ncol(indata)<100) 1.2 else 0.6, vertex.color=group.colors, main=metagene.filter.list[[i]]$n )
				legend( "bottomright", as.character(unique( group.labels )), cex=0.5, text.col=unique.group.colors, bg="white" )	
				box()

			plot( g, layout=layout, vertex.size=ifelse( ncol(indata)<250, 5, 3 ), vertex.label = rep("",ncol(indata)), vertex.label.cex=if(ncol(indata)<100) 1.2 else 0.6, vertex.color=group.colors, main=metagene.filter.list[[i]]$n )
				legend( "bottomright", as.character(unique( group.labels )), cex=0.5, text.col=unique.group.colors, bg="white" )	
				box()
			
			plot( g, layout=layout, vertex.size=ifelse( ncol(indata)<250, 5, 3 ), vertex.label = rep("",ncol(indata)), vertex.label.cex=if(ncol(indata)<100) 1.2 else 0.6, vertex.color=bootstrap.col, main=metagene.filter.list[[i]]$n )
				legend( "bottomright", as.character(unique( group.labels )), cex=0.5, text.col=unique.group.colors, bg="white" )	
				box()
			
			

# 			if( length(unique(group.labels)) > 1 )
# 			{	
# 				
# 				# collapse CN nodes
# 	
# 				adj.matrix = as.matrix( get.adjacency(g) )
# 				
# 				hypernode.strength = rep( 1, ncol(indata) )
# 				names(hypernode.strength) = colnames(indata)
# 				
# 				while( TRUE )
# 				{
# 	
# 	
# 					cliqs = maximal.cliques(g)
# 					for( ii in 1:length(cliqs) )
# 					{
# 						cliq.nodes = get.vertex.attribute(g,"name")[ cliqs[[ii]] ]
# 						max.group.label = names( which.max( table( group.labels[ cliq.nodes ] ) ) )
# 					
# 						cliqs[[ii]] = cliq.nodes[ which( group.labels[ cliq.nodes ] == max.group.label ) ]
# 						
# 						cliqs[[ii]] = names( which( hypernode.strength[cliqs[[ii]]] == 1 ) )
# 					}
# 						
# 					max.cliq = cliqs[[ which.max( sapply( cliqs, length ) ) ]]
# 				
# 					if( length(max.cliq) < 2 ) break
# 	
# 					
# 					adj.matrix[max.cliq[1],] = adj.matrix[max.cliq[1],] + apply(adj.matrix[max.cliq[-1],,drop=F],2,sum)
# 					adj.matrix[,max.cliq[1]] = adj.matrix[max.cliq[1],]
# 					diag(adj.matrix) = 0
# 					adj.matrix = adj.matrix[ -which( rownames(adj.matrix) %in% max.cliq[-1] ), -which( colnames(adj.matrix) %in% max.cliq[-1] ) ]
# 					
# 					hypernode.strength[ max.cliq[1] ] = hypernode.strength[ max.cliq[1] ] + sum( hypernode.strength[ max.cliq[-1] ] )
# 					hypernode.strength = hypernode.strength[ -which( names(hypernode.strength) %in% max.cliq[-1] ) ]
# 					
# 					g = graph.adjacency( adj.matrix, weighted=T,  mode="undirected" )
# 	
# 				}
# 	
# 	
# 				edge.width = get.edge.attribute(g, "weight")			
# 				
# 				plot( g, layout=layout[match( colnames(adj.matrix), colnames(indata) ),], vertex.size=5-(sqrt(1/max(table(group.labels)))*30)+sqrt(hypernode.strength/max(table(group.labels)))*30, vertex.label = hypernode.strength , vertex.label.cex=1.2, vertex.color=group.colors[rownames(adj.matrix)], vertex.label.color="white", edge.width=edge.width/max(edge.width)*10, main=metagene.filter.list[[i]]$n )			
# 					legend( "bottomright", as.character(unique( group.labels )), cex=0.5, text.col=unique.group.colors, bg="white" )	
# 					box()			
# 				
# 				plot( g, layout=layout.fruchterman.reingold, vertex.size=5-(sqrt(1/max(table(group.labels)))*30)+sqrt(hypernode.strength/max(table(group.labels)))*30, vertex.label = hypernode.strength , vertex.label.cex=1.2, vertex.color=group.colors[rownames(adj.matrix)], vertex.label.color="white", edge.width=edge.width/max(edge.width)*10, main=metagene.filter.list[[i]]$n )			
# 					legend( "bottomright", as.character(unique( group.labels )), cex=0.5, text.col=unique.group.colors, bg="white" )	
# 					box()			
#  		
# 			}
			
		
		}


		
		
		
		


		# PCM

		
		hcl = hclust( dist( cor(s) ) )		
		
		heatmap.wrap( x=cor( s ), Rowv=as.dendrogram(hcl), Colv=as.dendrogram(hcl), col=colramp(1000), scale="n", main=metagene.filter.list[[i]]$n, mar=c(8,8), ColSideColors=group.colors, RowSideColors=group.colors )
			par(new=T)
			plot(0,type="n", axes=F, xlab="", ylab="" )
			legend( "bottomright", as.character(unique( group.labels )), cex=0.5, text.col=unique.group.colors, bg="white" )	

		
		if( i <= 2 )
		{
			def.par <- par(no.readonly = TRUE)
			source("lib/f - 2nd lvl Module Correlation Map.r")
			par(def.par) 
		}
		
		
		heatmap.wrap( x=cor( s ), Rowv=NA, Colv=NA, col=colramp(1000), scale="n", main=metagene.filter.list[[i]]$n, mar=c(8,8), ColSideColors=group.colors, RowSideColors=group.colors )
			par(new=T)
			plot(0,type="n", axes=F, xlab="", ylab="" )
			legend( "bottomright", as.character(unique( group.labels )), cex=0.5, text.col=unique.group.colors, bg="white" )	



	}


	dev.off()





