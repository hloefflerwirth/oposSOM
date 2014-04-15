



	pdf( paste( files.name, "- Results/Rank Maps.pdf" ) , 29.7/2.54, 21/2.54 )


	par( mar=c(0,0,0,0), mfrow=c(1,1) )
	plot( 0, type="n", xlab="", ylab="", axes=T )
	text( 1, 0.3, "<FoldChange Rank> Maps", cex=2 )



	par( mfrow=c( 7, 12 ) )
	par( mar=c(0.3,0.9,4.5,0.9) )
	count.col = 0
	for( gl in 1:length( unique( group.labels ) ) )
	{
		plot( 0, type="n", axes=F, xlab="", ylab="", xlim=c(0,1) )
		if( length( unique( group.labels ) ) > 1 )
			mtext( unique( group.labels )[gl], side=3, line = 2, cex=1.5, at=0, font=3, adj=0, col=unique.group.colors[gl] )

		par(new=T)
		for( j in which( group.labels == unique( group.labels )[gl] ) )
		{
			mean.rank = tapply( indata[,j], som.nodes, mean ) 
			mean.rank.full = rep( NA, preferences$dim.som1^2 )
			names( mean.rank.full ) = as.character( 1:preferences$dim.som1^2 )
			mean.rank.full[ names(mean.rank) ] = mean.rank

			image( matrix( mean.rank.full, preferences$dim.som1, preferences$dim.som1 ), axes=F, col = colramp(1000) )
			title( paste(j,":",colnames(indata)[j]), line=1, cex.main=0.8 )
			box()
			

			count.col = count.col + 1 
		}
		if( count.col %% 12 != 0 )
		for( j in 1:( 12 - count.col %% 12 ) )
		{
			plot( 0, type="n", axes=F, xlab="", ylab="" )
			
			count.col = count.col + 1
		}
	}








	par( mar=c(0,0,0,0), mfrow=c(1,1) )
	plot( 0, type="n", xlab="", ylab="", axes=T )
	text( 1, 0.3, "<WAD Rank> Maps", cex=2 )



	par( mfrow=c( 7, 12 ) )
	par( mar=c(0.3,0.9,4.5,0.9) )
	WADv = c()
	count.col = 0
	for( gl in 1:length( unique( group.labels ) ) )
	{
		plot( 0, type="n", axes=F, xlab="", ylab="", xlim=c(0,1) )
		if( length( unique( group.labels ) ) > 1 )
			mtext( unique( group.labels )[gl], side=3, line = 2, cex=1.5, at=0, font=3, adj=0, col=unique.group.colors[gl] )

		par(new=T)
		for( j in which( group.labels == unique( group.labels )[gl] ) )
		{
			mean.rank = tapply( WAD.g.m[,j], som.nodes, mean ) 
			mean.rank.full = rep( NA, preferences$dim.som1^2 )
			names( mean.rank.full ) = as.character( 1:preferences$dim.som1^2 )
			mean.rank.full[ names(mean.rank) ] = mean.rank

			image( matrix( mean.rank.full, preferences$dim.som1, preferences$dim.som1 ), axes=F, col = colramp(1000) )
			title( paste(j,":",colnames(indata)[j]), line=1, cex.main=0.8 )
			box()

			count.col = count.col + 1 
		}
		if( count.col %% 12 != 0 )
		for( j in 1:( 12 - count.col %% 12 ) )
		{
			plot( 0, type="n", axes=F, xlab="", ylab="" )
			
			count.col = count.col + 1
		}
	}








	suppressWarnings (
	{


	par( mar=c(0,0,0,0), mfrow=c(1,1) )
	plot( 0, type="n", xlab="", ylab="", axes=T )
	text( 1, 0.3, "<Shrinkage-t Rank> Maps", cex=2 )



	par( mfrow=c( 7, 12 ) )
	par( mar=c(0.3,0.9,4.5,0.9) )
	WADv = c()
	count.col = 0
	for( gl in 1:length( unique( group.labels ) ) )
	{
		plot( 0, type="n", axes=F, xlab="", ylab="", xlim=c(0,1) )
		if( length( unique( group.labels ) ) > 1 )
			mtext( unique( group.labels )[gl], side=3, line = 2, cex=1.5, at=0, font=3, adj=0, col=unique.group.colors[gl] )

		par(new=T)
		for( j in which( group.labels == unique( group.labels )[gl] ) )
		{
			mean.rank = tapply( t.g.m[,j], som.nodes, mean ) 
			mean.rank.full = rep( NA, preferences$dim.som1^2 )
			names( mean.rank.full ) = as.character( 1:preferences$dim.som1^2 )
			mean.rank.full[ names(mean.rank) ] = mean.rank

			image( matrix( mean.rank.full, preferences$dim.som1, preferences$dim.som1 ), axes=F, col = colramp(1000) )
			title( paste(j,":",colnames(indata)[j]), line=1, cex.main=0.8 )
			box()

			count.col = count.col + 1 
		}
		if( count.col %% 12 != 0 )
		for( j in 1:( 12 - count.col %% 12 ) )
		{
			plot( 0, type="n", axes=F, xlab="", ylab="" )
			
			count.col = count.col + 1
		}
	}




	} )




	dev.off()










