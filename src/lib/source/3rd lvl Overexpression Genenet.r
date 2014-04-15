

require.cran( "GeneNet" )	





	spot.profiles = t( sapply( GS.infos.overexpression$spots, function(x) colMeans(indata[x$genes,,drop=F]) ) )

	o = capture.output( { 
		inferred.pcor = ggm.estimate.pcor(t(spot.profiles)); 
		suppressWarnings( { test.results = ggm.test.edges(inferred.pcor, plot=F, verbose=F) } )
	})

	
	sig.edges = which( test.results$prob > 0.9 )
	
	if( length(sig.edges) > 0 )
	{
		
		
		edge.list = test.results[sig.edges,]
		edge.list$node1 = rownames(spot.profiles)[ edge.list$node1 ]
		edge.list$regulates = rownames(spot.profiles)[ edge.list$node2 ]
		edge.list$edgetype = sign(edge.list$pcor)
		
		gene.list = by( edge.list[,c("regulates","prob","pcor","edgetype")], edge.list[,"node1"], c )
		
		
		in.only.genes = setdiff( edge.list$regulates, names(gene.list) )
		in.only.gene.list = lapply( in.only.genes, function(x) list( regulates=c(), pcor=c(), edgetype=c() ) )
		names(in.only.gene.list) = in.only.genes
		
		gene.list = c( gene.list, in.only.gene.list )
		
	
		
		
		
		adj.list.numeric = lapply( gene.list, function(x){ match( x$regulates, names(gene.list) )  } )				
		pcor.list = unlist( lapply( gene.list, function(x){ x$pcor } ) )
		
		g = graph.adjlist( adj.list.numeric )		
		g = set.vertex.attribute( g, "name", value=names(adj.list.numeric) )
		
		
		n.gene.spot = sapply( GS.infos.overexpression$spots, function(x) length(x$genes) )[ get.vertex.attribute(g,"name") ]
		vertex.size = 8*n.gene.spot/max(n.gene.spot)+2

		edge.width = 5*abs(pcor.list)/max(abs(pcor.list))+1
		edge.color = c()
		for(i in 1:length(gene.list) )
			if( length(gene.list[[i]]$regulates) > 0 )
			{
				vertex.edge.color = c("red","green")[ gene.list[[i]]$edgetype/2 + 1.5 ]
				undirected = sapply( gene.list[ gene.list[[i]]$regulates ], function(x) names(gene.list)[i] %in% x$regulates )
				vertex.edge.color[ undirected ] = "gray50"						
				edge.color = c( edge.color, vertex.edge.color )
			}
		

		
		
		
		n = paste( files.name, " - Results/3rd lvl Spot Analysis/Overexpression GeneNet", if(GS.infos.overexpression$filtered) " filtered", ".pdf", sep="" )	
		pdf( n, 29.7/2.54, 21/2.54 )
		

		par(mar=c(4.7,13.2,4.7,13.2))
		image( matrix(GS.infos.overexpression$overview.mask, preferences$dim.som1), col="gray90", axes=F )
		par(new=T, mar=c(1,1,1,1))
		plot( g, layout=t( sapply( GS.infos.overexpression$spots[get.vertex.attribute(g,"name")], function(x) x$position ) ), vertex.label=get.vertex.attribute(g,"name"), vertex.label.color="black" ,vertex.label.cex=1, vertex.color="grey", vertex.size=vertex.size, edge.color=edge.color, edge.width=edge.width, main="", xlim=c(-1.2,1.2), ylim=c(-1.2,1.2) )		
		rect(-1.1,-1.1,1.1,1.1)

		
		dev.off()
	
	}
