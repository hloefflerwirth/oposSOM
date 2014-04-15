


	supersom.custom = som( t(metadata), xdim=preferences$dim.som2, ydim=preferences$dim.som2 )



	if( preferences$dim.som2 == 20 )
	{
		supersom.20 = supersom.custom
	} else
	{
		supersom.20 = som( t(metadata), xdim=20, ydim=20 )
	}






	pdf( paste( files.name, "- Results/2nd lvl Metagene Analysis/2nd lvl SOM.pdf" ), 21/2.54, 21/2.54 )



	##### Plot Supersom #####

	par( mar=c(1,1,1,1) )




	xl = c( min( supersom.custom$visual[,"x"] )-1.2, max( supersom.custom$visual[,"x"] )+1.2 )
	yl = c( min( -supersom.custom$visual[,"y"] )-1.2, max( -supersom.custom$visual[,"y"] )+1.2 )
	plot( supersom.custom$visual[,"x"], -supersom.custom$visual[,"y"], type="n", axes=F, xlab="", ylab="", xlim=xl, ylim=yl, xaxs="i", yaxs="i"  )

	if( ncol(indata) < 100 )
	{
		legend( "bottomright",  paste( 1:length( colnames( indata ) ), ":", colnames( indata ) ), cex=0.5, text.col=group.colors, ncol=(ncol(indata)-1)%/%25+1, bg="white" )
	}
	if( length(unique( group.labels ) ) > 1 )
	{
		legend( "topright", as.character(unique( group.labels )), cex=0.5, text.col=unique.group.colors, bg="white" )
	}

	for( j in 1:nrow(supersom.custom$code.sum) )
	{

		which.samples = intersect( which( supersom.custom$visual[,"x"] == supersom.custom$code.sum[j,"x"] ),
						which( supersom.custom$visual[,"y"] == supersom.custom$code.sum[j,"y"] )	)

		if( !is.na( which.samples[1] ) )
		{

			which.samples = which.samples[ 1:min( 4, length(which.samples) ) ]


			x.seq = rep( c( -0.2, 0.2 ), 2 ) [ 1:length(which.samples) ]
			y.seq = c( rep( 0.2, 2 ), rep( -0.2, 2 ) ) [ 1:length(which.samples) ]


			points( supersom.custom$visual[ which.samples[1], "x" ]+x.seq, -supersom.custom$visual[ which.samples[1], "y" ]+y.seq,
				pch=16, col=group.colors[ which.samples ], cex=2.5 )

			points( supersom.custom$visual[ which.samples[1], "x" ]+x.seq, -supersom.custom$visual[ which.samples[1], "y" ]+y.seq,
				pch=1, col="gray20", cex=2.5, lwd=1 )

			text( supersom.custom$visual[ which.samples[1], "x" ]+x.seq, -supersom.custom$visual[ which.samples[1], "y" ]+y.seq,
				which.samples, col="gray20", cex=0.8 )

		}
	}

	box()




	plot( supersom.custom$visual[,"x"], -supersom.custom$visual[,"y"], type="n", axes=F, xlab="", ylab="", xlim=xl, ylim=yl, xaxs="i", yaxs="i" )

	if( length(unique( group.labels ) ) > 1 )
	{
		legend( "topright", as.character(unique( group.labels )), cex=0.5, text.col=unique.group.colors, bg="white" )
	}

	for( j in 1:nrow(supersom.custom$code.sum) )
	{

		which.samples = intersect( which( supersom.custom$visual[,"x"] == supersom.custom$code.sum[j,"x"] ),
															 which( supersom.custom$visual[,"y"] == supersom.custom$code.sum[j,"y"] )	)

		if( !is.na( which.samples[1] ) )
		{

			which.samples = which.samples[ 1:min( 4, length(which.samples) ) ]


			x.seq = rep( c( -0.2, 0.2 ), 2 ) [ 1:length(which.samples) ]
			y.seq = c( rep( 0.2, 2 ), rep( -0.2, 2 ) ) [ 1:length(which.samples) ]


			points( supersom.custom$visual[ which.samples[1], "x" ]+x.seq, -supersom.custom$visual[ which.samples[1], "y" ]+y.seq,
							pch=16, col=group.colors[ which.samples ], cex=2.5 )

			points( supersom.custom$visual[ which.samples[1], "x" ]+x.seq, -supersom.custom$visual[ which.samples[1], "y" ]+y.seq,
							pch=1, col="gray20", cex=2.5, lwd=1 )

			text( supersom.custom$visual[ which.samples[1], "x" ]+x.seq, -supersom.custom$visual[ which.samples[1], "y" ]+y.seq,
						colnames(indata)[which.samples], col="gray20", cex=0.8 )

		}
	}

	box()




	##### Polygone-2ndSOM #####

	if( length(unique( group.labels ) ) > 1 )
	{

		transparent.group.colors = sapply( unique.group.colors, function(x)
		{
			paste( substr( x, 1, 7 ) , "50", sep="" )
		} )
		names(transparent.group.colors) = unique(group.labels)



		plot( supersom.custom$visual[,"x"], -supersom.custom$visual[,"y"], type="n", axes=F, xlab="", ylab="", cex=4, col=group.colors, pch=16, xaxs="i", yaxs="i", xlim=xl, ylim=yl )


		for( i in 1:length(unique(group.labels)) )
		{

			group.member = which(group.labels==unique(group.labels)[i])
			group.centroid = colMeans( supersom.custom$visual[ group.member, 1:2 ])




			hull = chull(supersom.custom$visual[ group.member, 1 ], -supersom.custom$visual[ group.member, 2 ])
			polygon(supersom.custom$visual[ group.member[hull], 1 ], -supersom.custom$visual[ group.member[hull], 2 ], col=transparent.group.colors[i], lty=1, border=unique.group.colors[i])


			# 			m = matrix( group.metadata[,i], preferences$dim.som1, preferences$dim.som1 )
			# 			if( max(m) - min(m) != 0 )
			# 				m = ( m - min(m) ) / ( max(m) - min(m) )*999
			# 			m = cbind( apply( m, 1, function(x){x} ) )[ nrow(m):1, ]
			#
			# 			x <- pixmapIndexed( m , col = colramp(1000) )
			# 			addlogo( x, group.centroid[1]+preferences$dim.som2*c(-0.04,0.04), -group.centroid[2]+preferences$dim.som2*c(-0.04,0.04) )
			# 			rect( group.centroid[1]-preferences$dim.som2*0.04, -group.centroid[2]-preferences$dim.som2*0.04, group.centroid[1]+preferences$dim.som2*0.04, -group.centroid[2]+preferences$dim.som2*0.04 )
			#
			#
			# 			text( group.centroid[1], -group.centroid[2]+preferences$dim.som2*0.08, unique(group.labels)[i], cex=1.6, col=unique.group.colors[i] )
			#


		}



		legend( "topright", as.character(unique( group.labels )), cex=0.5, text.col=unique.group.colors, bg="white" )

		box()


	}



	##### Plot SmoothSupersom ######

	source("R/smooth_scatterplot.r")

	for( i in 1:length(unique(group.labels)) )
	{
		col = colorRampPalette(c("white", unique.group.colors[i] ))
		coords = cbind( supersom.custom$visual[,"x"], -supersom.custom$visual[,"y"] )
		coords = coords[ which(group.labels==unique(group.labels)[i]) , ,drop=F]
		if( nrow(coords) == 1 ) coords = rbind(coords,coords)

		smoothScatter( coords, main="", xlim=xl, ylim=yl, axes=F, xlab="", ylab="", nrpoints=0, colramp=col, bandwidth=2, nbin=128, transformation = function(x) x^0.25, xaxs="i", yaxs="i" )

		if( i < length(unique(group.labels)) ) par(new=T)
	}
	par(new=T)
	plot( supersom.custom$visual[,"x"], -supersom.custom$visual[,"y"], pch=16, col=group.colors, axes=F, xlab="",ylab="", xlim=xl,ylim=yl, xaxs="i", yaxs="i" )

	if( length(unique( group.labels ) ) > 1 )
	{
		legend( "topright", as.character(unique( group.labels )), cex=0.5, text.col=unique.group.colors, bg="white" )
	}





	##### Plot Supersom with real expression profiles ######



	par( mar=c(1,1,1,1) )



	xl = c( min( supersom.20$visual[,"x"] )-1, max( supersom.20$visual[,"x"] )+1 )
	yl = c( min( -supersom.20$visual[,"y"] )-1, max( -supersom.20$visual[,"y"] )+1 )
	plot( supersom.20$visual[,"x"], -supersom.20$visual[,"y"], type="p", axes=F, xlab="", ylab="", xlim=xl, ylim=yl, xaxs="i", yaxs="i" )

	if( ncol(indata) < 100 )
	{
		legend( "bottomright",  paste( 1:length( colnames( indata ) ), ":", colnames( indata ) ), cex=0.5, text.col=group.colors, ncol=(ncol(indata)-1)%/%25+1, bg="white" )
	}
	if( length(unique( group.labels ) ) > 1 )
	{
		legend( "topright", as.character(unique( group.labels )), cex=0.5, text.col=unique.group.colors, bg="white" )
	}

	for( j in 1:nrow(supersom.20$code.sum) )
	{

		which.samples = intersect( which( supersom.20$visual[,"x"] == supersom.20$code.sum[j,"x"] ),
						which( supersom.20$visual[,"y"] == supersom.20$code.sum[j,"y"] )	)

		if( !is.na( which.samples[1] ) )
		{

			m = matrix( metadata[, which.samples[1] ], preferences$dim.som1, preferences$dim.som1 )
			if( max(m) - min(m) != 0 )
				m = ( m - min(m) ) / ( max(m) - min(m) )*999
			m = cbind( apply( m, 1, function(x){x} ) )[ nrow(m):1, ]

			x <- pixmapIndexed( m , col = colramp(1000) )
			addlogo( x, supersom.20$visual[ which.samples[1], "x" ]+c(-0.45,0.455), -supersom.20$visual[ which.samples[1], "y" ]+c(-0.45,0.45) )


			which.samples = which.samples[ 1:min( 4, length(which.samples) ) ]


			x.seq = rep( c( -0.2, 0.2 ), 2 ) [ 1:length(which.samples) ]
			y.seq = c( rep( 0.2, 2 ), rep( -0.2, 2 ) ) [ 1:length(which.samples) ]


			points( supersom.20$visual[ which.samples[1], "x" ]+x.seq, -supersom.20$visual[ which.samples[1], "y" ]+y.seq,
				pch=16, col=group.colors[ which.samples ], cex=2.5 )

			points( supersom.20$visual[ which.samples[1], "x" ]+x.seq, -supersom.20$visual[ which.samples[1], "y" ]+y.seq,
				pch=1, col="gray20", cex=2.5, lwd=1 )

			text( supersom.20$visual[ which.samples[1], "x" ]+x.seq, -supersom.20$visual[ which.samples[1], "y" ]+y.seq,
				which.samples, col="gray20", cex=0.8 )

		}
	}

	box()














	dev.off()





