

	verbose = T
	if( output.paths["LPE"] == "" ) verbose = F

	
	if(verbose)
		write.table( preferences$error.model, paste( output.paths["LPE"],"/1_error_model.txt", sep="" ), row.names=F, col.names=F )
	
	

	cat("Single Genes\n|" )
	for( i in 1:50 ) cat(" ")
	cat("|\n|")
	flush.console()





	mean.e.g = rowMeans( indata.original )
	e.g.m = matrix( NA, nrow(indata), ncol(indata), dimnames=list( rownames( indata ), colnames( indata ) ) )
	delta.e.g.m = matrix( NA, nrow(indata), ncol(indata), dimnames=list( rownames( indata ), colnames( indata ) ) )

	sd.g.m = matrix( NA, nrow(indata), ncol(indata), dimnames=list( rownames( indata ), colnames( indata ) ) )
	LPE.g.m = matrix( NA, nrow(indata), ncol(indata), dimnames=list( rownames( indata ), colnames( indata ) ) )

	WAD.g.m = matrix( NA, nrow(indata), ncol(indata), dimnames=list( rownames( indata ), colnames( indata ) ) )

	t.g.m = matrix( NA, nrow(indata), ncol(indata), dimnames=list( rownames( indata ), colnames( indata ) ) )
	p.g.m = matrix( NA, nrow(indata), ncol(indata), dimnames=list( rownames( indata ), colnames( indata ) ) )


	n.0.m = rep( NA, ncol(indata) )
	names(n.0.m) = colnames(indata)
	perc.DE.m = rep( NA, ncol(indata) )
	names(perc.DE.m) = colnames(indata)

	fdr.g.m = matrix( NA, nrow(indata), ncol(indata), dimnames=list( rownames( indata ), colnames( indata ) ) )
	Fdr.g.m = matrix( NA, nrow(indata), ncol(indata), dimnames=list( rownames( indata ), colnames( indata ) ) )

	mean.LPE2 = rep( NA, ncol(indata) )
	names(mean.LPE2) = colnames(indata)




	for( m in 1:ncol(indata)  ) 
	{
		if( preferences$error.model == "replicates" )	
		{
			e.r.g.m = as.matrix( indata.original[ , which( colnames(indata.original) == unique( colnames(indata) )[m] ) ,drop=F ] )
		} else
		{
			e.r.g.m = as.matrix( indata.original[ , m ] )
		} 



		e.g.m[ , m ] = rowMeans( e.r.g.m ) 

		delta.e.r.g.m = e.r.g.m - mean.e.g
		delta.e.g.m[,m] = rowMeans( delta.e.r.g.m )

		w.g.m = ( delta.e.g.m[,m] - min(delta.e.g.m[,m]) ) / ( max(delta.e.g.m[,m]) - min(delta.e.g.m[,m]) )
		WAD.g.m[ , m ] = w.g.m * delta.e.g.m[,m]
	}		
	





	if( preferences$error.model == "replicates" )			##################################
	{
		R.m = c()

		for( m in 1:ncol(indata) ) 
		{
			R.m[m] = length( which( colnames( indata.original ) == colnames( indata )[m] ) )
		}
		names( R.m ) = colnames( indata )


		sd.shrink.g.m = matrix( NA, nrow(indata), ncol(indata), dimnames=list( rownames( indata ), colnames( indata ) ) )

		for( m in 1:ncol(indata)  ) 
		{
			if( R.m[ m ] > 1 )
			{
				e.r.g.m = as.matrix( indata.original[ , which( colnames(indata.original) == unique( colnames(indata) )[m] ) ,drop=F ] )

				o = order( e.g.m[,m] )

				sd.g.m[ , m ] = sapply( c(1:nrow(indata)), function(x)
				{
					sqrt( sum( ( e.r.g.m[x,]-e.g.m[x,m])^2 )/R.m[m] ) 
				} )


				LPE.g.m[ o , m ] = Get.Running.Average( sd.g.m[o,m], min( 200, round( nrow(indata)*0.02 ) ) )


#				if( exists( "chip.LPE.monotony" ) )
				{
#					if( chip.LPE.monotony == "falling" )
					{
						SD2 = LPE.g.m[ o , m ]
						for( j in (length(SD2)-1):1 )
						{		
							if( SD2[ j ] < SD2[ j+1 ] ) SD2[ j ] = SD2[ j+1 ]
						}
						LPE.g.m[ o , m ] = SD2
					}
				}
	
				lambda = 0.5
				sd.shrink.g.m[ , m ] = sqrt( lambda * sd.g.m[ , m ]^2 + ( 1 - lambda ) * LPE.g.m[ , m ]^2 )
			}

		out.intervals = round( seq( 1, ncol(indata), length.out=26 ) )[-1]
		cat( paste( rep("#",length( which( out.intervals == m ) ) ), collapse="" ) );	flush.console()
		}


		no.sd.samples = colnames(indata)[ which( is.na(apply( sd.shrink.g.m, 2, mean )) )  ]
		if( length(no.sd.samples) > 0 )
		{
			for( i in 1:nrow(indata) )
			{

				gene.expression.order = colnames( indata )[ order( indata[i,] ) ]

				for( m in no.sd.samples )
				{
					position = which( gene.expression.order == m )

					gene.expression.order.help = gene.expression.order
					gene.expression.order.help[ gene.expression.order.help %in% no.sd.samples ] = NA
					gene.expression.order.help[ position ] = m
					gene.expression.order.help = na.omit( gene.expression.order.help )
	
					position = which( gene.expression.order.help == m )
			
					window = position + c( -2, -1, 1, 2 )
					window[ which( window < 1 ) ] = NA
					window[ which( window > length(gene.expression.order.help) ) ] = NA
					window = na.omit( window )
			
					sd.g.m[ i, m ] = mean( sd.g.m[ i, gene.expression.order.help[ window ] ] )
				}
			}

			for( m in no.sd.samples )
			{
				o = order( e.g.m[,m] )

				LPE.g.m[ o , m ] = Get.Running.Average( sd.g.m[o,m], min( 200, round( nrow(indata)*0.02 ) ) )

#				if( exists( "chip.LPE.monotony" ) )
				{
#					if( chip.LPE.monotony == "falling" )
					{
						SD2 = LPE.g.m[ o , m ]
						for( j in (length(SD2)-1):1 )
						{		
							if( SD2[ j ] < SD2[ j+1 ] ) SD2[ j ] = SD2[ j+1 ]
						}
						LPE.g.m[ o , m ] = SD2
					}
				}
	
				lambda = 0.5
				sd.shrink.g.m[ , m ] = sqrt( lambda * sd.g.m[ , m ]^2 + ( 1 - lambda ) * LPE.g.m[ , m ]^2 )
			}
		} # END no.sd.samples



		if( any( sd.shrink.g.m == 0 ) )
		{
			sd.shrink.g.m[which(sd.shrink.g.m==0)] = max(sd.shrink.g.m[-which(sd.shrink.g.m==0)])

			cat( "\n!!Data contains constant genes!!\n" ); flush.console()			
			cat( "!!Constant genes set to maximum variance!!\n\n" ); flush.console()			
		}


		for( m in 1:ncol(indata) ) 
		{
			t.g.m[,m] = sqrt(R.m[m]) * ( e.g.m[,m ] - mean.e.g ) / sd.shrink.g.m[,m]

		out.intervals = round( seq( 1, ncol(indata), length.out=5+1 ) )[-1]
		cat( paste( rep("#",length( which( out.intervals == m ) ) ), collapse="" ) );	flush.console()
		}


	} # END error.model "replicates"






	if( preferences$error.model == "all.samples" )			##################################
	{
		for( m in 1:ncol(indata)  ) 
		{
			o = order( apply( indata.original, 1, mean ) )
		
			sd = apply( indata.original, 1, sd )

			LPE.g.m[ o , m ] = Get.Running.Average( sd[o] , min( 200, round( nrow(indata)*0.02 ) ) )
			LPE.g.m[ which( is.nan( LPE.g.m[ , m ] ) ), m ] = 0.0000000001
			LPE.g.m[ which( LPE.g.m[ , m ] == 0 ), m ] = 0.0000000001


#			if( exists( "chip.LPE.monotony" ) )
			{
#				if( chip.LPE.monotony == "falling" )
				{
					SD2 = LPE.g.m[ o , m ]
					for( j in (length(SD2)-1):1 )
					{		
						if( SD2[ j ] < SD2[ j+1 ] ) SD2[ j ] = SD2[ j+1 ]
					}
					LPE.g.m[ o , m ] = SD2
					rm(SD2)
				}
			}

		out.intervals = round( seq( 1, ncol(indata), length.out=26 ) )[-1]
		cat( paste( rep("#",length( which( out.intervals == m ) ) ), collapse="" ) );	flush.console()
		}


		sd.g.m = LPE.g.m


		for( m in 1:ncol(indata) ) 
		{
				t.g.m[,m] = sqrt(ncol(indata)) * ( e.g.m[,m ] - mean.e.g ) / sd.g.m[,m]

		out.intervals = round( seq( 1, ncol(indata), length.out=5+1 ) )[-1]
		cat( paste( rep("#",length( which( out.intervals == m ) ) ), collapse="" ) );	flush.console()
		}		

	} # END error.model "all.samples"




	if( preferences$error.model == "groups" )			##################################
	{
		R.m = c()

		for( m in 1:ncol(indata) ) 
		{
			R.m[m] = length( which( group.labels == group.labels[m] ) )
		}
		names( R.m ) = colnames( indata )


		sd.shrink.g.m = matrix( NA, nrow(indata), ncol(indata), dimnames=list( rownames( indata ), colnames( indata ) ) )

		for( m in 1:ncol(indata)  ) 
		{
			if( R.m[ m ] > 1 )
			{
				e.r.g.m = as.matrix( indata.original[ , which( group.labels == group.labels[m] ) ,drop=F ] )

				o = order( e.g.m[,m] )

				sd.g.m[ , m ] = sapply( c(1:nrow(indata)), function(x)
				{
					sqrt( sum( ( e.r.g.m[x,]-e.g.m[x,m])^2 )/R.m[m] ) 
				} )


				LPE.g.m[ o , m ] = Get.Running.Average( sd.g.m[o,m], min( 200, round( nrow(indata)*0.02 ) ) )


#				if( exists( "chip.LPE.monotony" ) )
				{
#					if( chip.LPE.monotony == "falling" )
					{
						SD2 = LPE.g.m[ o , m ]
						for( j in (length(SD2)-1):1 )
						{		
							if( SD2[ j ] < SD2[ j+1 ] ) SD2[ j ] = SD2[ j+1 ]
						}
						LPE.g.m[ o , m ] = SD2
					}
				}
	
				lambda = 0.5
				sd.shrink.g.m[ , m ] = sqrt( lambda * sd.g.m[ , m ]^2 + ( 1 - lambda ) * LPE.g.m[ , m ]^2 )
			}

		out.intervals = round( seq( 1, ncol(indata), length.out=26 ) )[-1]
		cat( paste( rep("#",length( which( out.intervals == m ) ) ), collapse="" ) );	flush.console()
		}


		no.sd.samples = colnames(indata)[ which( is.na(apply( sd.shrink.g.m, 2, mean )) )  ]
		if( length(no.sd.samples) > 0 )
		{
			for( i in 1:nrow(indata) )
			{

				gene.expression.order = colnames( indata )[ order( indata[i,] ) ]

				for( m in no.sd.samples )
				{
					position = which( gene.expression.order == m )

					gene.expression.order.help = gene.expression.order
					gene.expression.order.help[ gene.expression.order.help %in% no.sd.samples ] = NA
					gene.expression.order.help[ position ] = m
					gene.expression.order.help = na.omit( gene.expression.order.help )
	
					position = which( gene.expression.order.help == m )
			
					window = position + c( -2, -1, 1, 2 )
					window[ which( window < 1 ) ] = NA
					window[ which( window > length(gene.expression.order.help) ) ] = NA
					window = na.omit( window )
			
					sd.g.m[ i, m ] = mean( sd.g.m[ i, gene.expression.order.help[ window ] ] )
				}
			}

			for( m in no.sd.samples )
			{
				o = order( e.g.m[,m] )

				LPE.g.m[ o , m ] = Get.Running.Average( sd.g.m[o,m], min( 200, round( nrow(indata)*0.02 ) ) )

#				if( exists( "chip.LPE.monotony" ) )
				{
#					if( chip.LPE.monotony == "falling" )
					{
						SD2 = LPE.g.m[ o , m ]
						for( j in (length(SD2)-1):1 )
						{		
							if( SD2[ j ] < SD2[ j+1 ] ) SD2[ j ] = SD2[ j+1 ]
						}
						LPE.g.m[ o , m ] = SD2
					}
				}
	
				lambda = 0.5
				sd.shrink.g.m[ , m ] = sqrt( lambda * sd.g.m[ , m ]^2 + ( 1 - lambda ) * LPE.g.m[ , m ]^2 )
			}
		} # END no.sd.samples




		for( m in 1:ncol(indata) ) 
		{
				t.g.m[,m] = sqrt(R.m[m]) * ( e.g.m[,m ] - mean.e.g ) / sd.shrink.g.m[,m]

		out.intervals = round( seq( 1, ncol(indata), length.out=5+1 ) )[-1]
		cat( paste( rep("#",length( which( out.intervals == m ) ) ), collapse="" ) );	flush.console()
		}


	} # END error.model "groups"




















	### calculate significances ###


#	ylim.max = 0
	for( m in 1:ncol(indata) ) 
	{
#		p.g.m[,m] = 2 * pt( -abs( t.g.m[,m] ), 2 )
#		fdrtool.result = fdrtool( p.g.m[,m], statistic="pvalue", verbose=F, plot=F )
		
		
		suppressWarnings({ try.res = try( { fdrtool.result = fdrtool( t.g.m[,m], verbose=F, plot=F ) }, silent=T ) })
		
		if( class(try.res) != "try-error" )
		{								
			p.g.m[,m] = fdrtool.result$pval
			fdr.g.m[,m] = fdrtool.result$lfdr
			Fdr.g.m[,m] = fdrtool.result$qval
	
			n.0.m[ m ] = fdrtool.result$param[1,"eta0"]
			perc.DE.m[ m ] = 1 - n.0.m[ m ]
			
		} else # happens for eg phenotype data
		{					
			p.g.m[,m] = order( indata[,m] ) / nrow(indata)
			fdr.g.m[,m] = p.g.m[,m]
			Fdr.g.m[,m] = p.g.m[,m]
			
			n.0.m[ m ] = 0.5
			perc.DE.m[ m ] = 1 - n.0.m[ m ]
		}
			

#		h = hist( p.g.m[,m], bre=20, plot=F )
#		y.max = max( h$density)
#		if( y.max > ylim.max ) ylim.max = y.max

		mean.LPE2[ m ] = sqrt( Get.Area( e.g.m[ , m ], LPE.g.m[ , m ]^2  ) / ( max(e.g.m[ , m ]) - min(e.g.m[ , m ])) )


	out.intervals = round( seq( 1, ncol(indata), length.out=20+1 ) )[-1]
	cat( paste( rep("#",length( which( out.intervals == m ) ) ), collapse="" ) );	flush.console()
	}





	### LPE plots ###

	if( verbose )
	{
		
		if( preferences$error.model == "replicates" )
		{
			for( m in 1:ncol(indata) ) 
			{
				bmp( paste( output.paths["LPE"],"/", colnames(indata)[m], ".bmp", sep="" ), 600, 600 )
					par( mar=c( 5,6,4,5 ) )
					plot( sd.g.m[ , m ] ~ e.g.m[ , m ], xlab="e", ylab="", main=colnames(indata)[m], xlim=c(min(e.g.m,na.rm=T),max(e.g.m,na.rm=T)), ylim=c(0,max(sd.g.m,na.rm=T)), las=1, cex.main=2.5, cex.lab=2, cex.axis=2 )
					mtext( expression(sigma), side=2, line=4, cex= 2, las=2 )
	
					a = round( mean.LPE2[colnames(indata)[m]] ,2  )
					expr = paste(  "paste(\"<\", sigma[LPE] ,\"> = ",a," \", sep=\"\" )", sep="" )
					expr = paste("expression(", expr, ")", sep = "")
					mtext( eval(parse(text = expr)), line=-1.7, cex=2 )
	
					points( LPE.g.m[ , m ] ~ e.g.m[ , m ], col="green", pch=16 )
	
				dev.off()
			}	
	
		}else
		if( preferences$error.model == "all.samples" )
		{
			bmp( paste( output.paths["LPE"],"/all samples.bmp", sep="" ), 600, 600 )
				par( mar=c( 5,6,4,5 ) )
	
				plot( apply( indata.original, 1, sd ) ~ apply( indata.original, 1, mean ), xlab=expression(e[g]), ylab="", main="LPE", las=1, cex.main=2.5, cex.lab=2, cex.axis=2 )
					mtext( expression(sigma[g]), side=2, line=4, cex= 2, las=2 )
					points( LPE.g.m[ , 1 ] ~ apply( indata.original, 1, mean ), col="green", pch=16 )
		
			dev.off()
	
		}else
		if( preferences$error.model == "groups" )
		{
			for( m in 1:length(unique(group.labels)) ) 
			{
				plot.sample = which( group.labels == unique(group.labels)[m] )[1]
	
				bmp( paste( output.paths["LPE"],"/", unique(group.labels)[m], ".bmp", sep="" ), 600, 600 )
					par( mar=c( 5,6,4,5 ) )
					plot( sd.g.m[ , plot.sample ] ~ e.g.m[ , plot.sample ], xlab="e", ylab="", main=unique(group.labels)[m], xlim=c(min(e.g.m,na.rm=T),max(e.g.m,na.rm=T)), ylim=c(0,max(sd.g.m,na.rm=T)), las=1, cex.main=2.5, cex.lab=2, cex.axis=2 )
						mtext( expression(sigma), side=2, line=4, cex= 2, las=2 )
						points( LPE.g.m[ , plot.sample ] ~ e.g.m[ , plot.sample ], col="green", pch=16 )
	
				dev.off()
			}	
	
		}

	}










	
	
	
	





	cat("|\nMetagenes\n|" )
	for( i in 1:50 ) cat(" ")
	cat("|\n|")
	flush.console()




	t.m = matrix( NA, preferences$dim.som1^2, ncol(indata), dimnames=list( 1:preferences$dim.som1^2, colnames( indata ) ) )
	p.m = matrix( NA, preferences$dim.som1^2, ncol(indata), dimnames=list( 1:preferences$dim.som1^2, colnames( indata ) ) )

	fdr.m = matrix( NA, preferences$dim.som1^2, ncol(indata), dimnames=list( 1:preferences$dim.som1^2, colnames( indata ) ) )


	
	
	t.m.help = do.call( rbind, by( t.g.m, som.nodes, colMeans ) )	
	t.m[rownames(t.m.help),] = t.m.help
	
	cat( paste( rep("#",20 ), collapse="" ) );	flush.console()
	

	for( m in 1:ncol(indata)  ) 
	{
		
# 			for( j in 1:preferences$dim.som1^2 )
# 			{
# 				genes = names( which( som.nodes == j ) )
# 
# 				t.m[ j, m ] = mean( t.g.m[ genes, m ] )
#						
#					p.m[ j, m ] = 10 ^ mean( log10( p.g.m[ genes, m ] ) )
#					fdr.m[ j, m ] = mean( fdr.g.m[ genes, m ] )
#				}
		

		suppressWarnings({ try.res = try( { fdrtool.result = fdrtool( as.vector(na.omit(t.m[,m])), verbose=F, plot=F ) }, silent=T ) })
		
		if( class(try.res) != "try-error" )
		{								
			p.m[ which( !is.na(t.m[,m]) ), m ] = fdrtool.result$pval
			fdr.m[ which( !is.na(t.m[,m]) ), m ] = fdrtool.result$lfdr
			
		} else # happens for eg phenotype data
		{
			p.m[ which( !is.na(t.m[,m]) ), m ] = t.m[ which( !is.na(t.m[,m]) ), m ] / max( t.m[,m], na.rm=T )
			fdr.m[ which( !is.na(t.m[,m]) ), m ] = p.m[ which( !is.na(t.m[,m]) ), m ]
		}

				
			

	out.intervals = round( seq( 1, length( unique( colnames(indata) ) ), length.out=30+1 ) )[-1]
	cat( paste( rep("#",length( which( out.intervals == m ) ) ), collapse="" ) );	flush.console()
	}



	cat("|\n\n")
	flush.console()










	if( verbose )
	{
	
		pdf( paste( output.paths["LPE"],"/Sigma_LPE.pdf", sep="" ), 8, 8 )
	
		par( mar=c( 10,6,4,5 ) )
	
	
		mean.LPE.string = expression( paste("<", sigma[LPE], ">", sep="" ) )
	
	
		barplot( mean.LPE2, col=group.colors, main=mean.LPE.string, las=2, cex.main=2.5, cex.lab=2, cex.axis=2 )
			box()
	
	
			
		if( length( unique(group.labels) ) > 1 )
		{
	
			mean.LPE2.boxes = by( mean.LPE2, group.labels[ names( mean.LPE2 ) ], c )[ unique( group.labels ) ]
	
			boxplot( mean.LPE2.boxes, col=unique(group.colors), las=2, main=mean.LPE.string, cex.main=2.5, cex.axis=2, xaxt="n" )
				axis( 1, 1:length(mean.LPE2.boxes), names(mean.LPE2.boxes), las=2 )
	
		}
	
		dev.off()

	}



