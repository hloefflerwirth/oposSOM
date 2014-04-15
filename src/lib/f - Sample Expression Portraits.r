



	### Portraits


	pdf( paste( files.name, "- Results/Expression Portraits.pdf" ) , 29.7/2.54, 21/2.54 )


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
			image( matrix( metadata[,j], preferences$dim.som1, preferences$dim.som1 ), axes=F, col = colramp(1000) )
				title( paste(j,":",colnames(indata)[j]), line=1, cex.main=0.8 )
				if( length(unique( group.labels ) ) > 1 )	title( paste(group.bootstrap.score[j],"%",sep=""), line=0.2, cex.main=0.6, col.main=unique.group.colors[gl] )			
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
			md = matrix(metadata[,j],preferences$dim.som1,preferences$dim.som1)
			nrz = nrow(md)
				
			zfacet <- md[-1, -1] + md[-1, -nrz] + md[-nrz, -1] + md[-nrz, -nrz]  
			facetcol <- cut(zfacet, 1000)
		
			persp( md, axes=F, border=NA, expand=0.5, col=colramp(1000)[facetcol], phi=45, theta=-5, xlab="",ylab="",zlab="", box=T )
				title( paste(j,":",colnames(indata)[j]), line=1, cex.main=0.8 )

			count.col = count.col + 1 
		}
		if( count.col %% 12 != 0 )
		for( j in 1:( 12 - count.col %% 12 ) )
		{
			plot( 0, type="n", axes=F, xlab="", ylab="" )
			
			count.col = count.col + 1
		}
	}
	

	dev.off()




	






	pdf( paste( files.name, "- Results/Expression Portraits alternative.pdf" ) , 29.7/2.54, 21/2.54 )



	
	
	
	bleached.metadata = metadata
	
	for( i in 1:ncol(metadata) )
	{
		pos.metagenes = which( metadata[,i] >= 0 )
		neg.metagenes = which( metadata[,i] < 0 )
	
		bleached.metadata[pos.metagenes,i] = bleached.metadata[pos.metagenes,i,drop=F] - pmin(   bleached.metadata[pos.metagenes,i,drop=F]   ,   apply( metadata[pos.metagenes,-i,drop=F], 1, max )   )
		bleached.metadata[neg.metagenes,i] = bleached.metadata[neg.metagenes,i,drop=F] - pmax(   bleached.metadata[neg.metagenes,i,drop=F]   ,   apply( metadata[neg.metagenes,-i,drop=F], 1, min )   )
	}	




	par( mar=c(0,0,0,0), mfrow=c(1,1) )
	plot( 0, type="n", xlab="", ylab="", axes=T )
	text( 1, 0.3, "Sample Specific SOM Portraits", cex=2 )

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
			image( matrix( bleached.metadata[,j], preferences$dim.som1, preferences$dim.som1 ), axes=F, col = colramp(1000), zlim=range(bleached.metadata) )
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
	text( 1, 0.3, "Absolute Metagene Portraits", cex=2 )

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
			image( matrix( metadata[,j], preferences$dim.som1, preferences$dim.som1 ), axes=F, col = colramp(1000), zlim=c(min(metadata),max(metadata)) )
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
	text( 1, 0.3, "Significance Portraits", cex=2 )



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
			image( matrix( -p.m[,j], preferences$dim.som1, preferences$dim.som1 ), axes=F, col = colorRampPalette( c( "blue4","blue4","blue3","blue3","blue2","blue2","blue1","lightblue","darkgreen","#008B00","green3","green","yellow","gold","orange","red","darkred" ) )(1000), zlim=c(-1,0) )			
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
	text( 1, 0.3, "WAD Portraits", cex=2 )



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
			image( matrix( WAD.metadata[,j], preferences$dim.som1, preferences$dim.som1 ), axes=F, col = colramp(1000), cex.main=0.6 )
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
	text( 1, 0.3, "loglog FC Portraits", cex=2 )



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
			image( matrix( loglog.metadata[,j], preferences$dim.som1, preferences$dim.som1 ), axes=F, col = colramp(1000) )
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
	text( 1, 0.3, "Quantile Scale Portraits", cex=2 )

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
			md = matrix( metadata[,j], preferences$dim.som1, preferences$dim.som1 )								
			image( matrix( 1:(preferences$dim.som1)^2, preferences$dim.som1, preferences$dim.som1 ), axes=F, col = colramp(1000)[cut(md, quantile(md, seq(0,1, len = 1001)),include.lowest = TRUE)] )
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
	text( 1, 0.3, "Shrinkage-t Portraits ( log(t) )", cex=2 )



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
			t = log10( abs( t.m[,j] ) )
			t = t - min( t, na.rm=T )
			t = t * sign( t.m[,j] )
			image( matrix( t, preferences$dim.som1, preferences$dim.som1 ), axes=F, col = colramp(1000) )
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










