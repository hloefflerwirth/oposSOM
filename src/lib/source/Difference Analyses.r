
get.SD.estimate = function( data, samples, lambda=0.5 )
{
	
	if( length(samples) > 1 )
	{
		sd.g.m = apply( data[,samples], 1, sd )	
		LPE.g.m = rep( NA, nrow(data) ) 
		names(LPE.g.m) = rownames(data)
		
		o = order( rowMeans( data[,samples] ) )
		LPE.g.m[ o ] = Get.Running.Average( sd.g.m[o], min( 200, round( nrow(data)*0.02 ) ) )
		
		SD2 = LPE.g.m[ o ]
		for( j in (length(SD2)-1):1 )
		{		
			if( SD2[ j ] < SD2[ j+1 ] ) SD2[ j ] = SD2[ j+1 ]
		}
		LPE.g.m[ o ] = SD2
		
		sd.estimate = sqrt( lambda * sd.g.m^2 + ( 1 - lambda ) * LPE.g.m^2 )						
		
	} else
	{		
		sd = apply( data, 1, sd )
		LPE.g.m = rep( NA, nrow(data) ) 
		names(LPE.g.m) = rownames(data)
		
		o = order( rowMeans( data ) )				
		LPE.g.m[ o ] = Get.Running.Average( sd[o] , min( 200, round( nrow(data)*0.02 ) ) )
		LPE.g.m[ which( is.nan( LPE.g.m ) ) ] = 0.0000000001
		LPE.g.m[ which( LPE.g.m == 0 ) ] = 0.0000000001
				
		SD2 = LPE.g.m[ o ]
		for( j in (length(SD2)-1):1 )
		{		
			if( SD2[ j ] < SD2[ j+1 ] ) SD2[ j ] = SD2[ j+1 ]
		}
		LPE.g.m[ o ] = SD2
		
		sd.estimate = LPE.g.m
	}	
	
	return( sd.estimate )
}





perform.differences.analyses = function()
{	

	cat( "Differences Analyses\n\n" ); flush.console()

	dir.create( paste( files.name, "- Results/Summary Sheets - Differences" ), showWarnings=F )		
	dir.create( paste( files.name, "- Results/Summary Sheets - Differences/CSV Sheets" ), showWarnings=F )
	
	
	
	WAD.g.m = matrix( NA, nrow(indata), length(preferences$differences.list), dimnames=list( rownames(indata), names(preferences$differences.list) ) )
	t.g.m = matrix( NA, nrow(indata), length(preferences$differences.list), dimnames=list( rownames(indata), names(preferences$differences.list) ) )	
	p.g.m = matrix( NA, nrow(indata), length(preferences$differences.list), dimnames=list( rownames(indata), names(preferences$differences.list) ) )
 	fdr.g.m = matrix( NA, nrow(indata), length(preferences$differences.list), dimnames=list( rownames(indata), names(preferences$differences.list) ) )
 	Fdr.g.m = matrix( NA, nrow(indata), length(preferences$differences.list), dimnames=list( rownames(indata), names(preferences$differences.list) ) )
 	n.0.m = rep( NA, length(preferences$differences.list) )
	names(n.0.m) = names(preferences$differences.list)
 	perc.DE.m = rep( NA, length(preferences$differences.list) )
 	names(perc.DE.m) = names(preferences$differences.list)		
	
	for( d in 1:length(preferences$differences.list) )
	{
		
		samples.indata = list( preferences$differences.list[[d]][[1]], preferences$differences.list[[d]][[2]] )
		samples.indata.original = list(  which( colnames(indata.original) %in% colnames(indata)[ samples.indata[[1]] ] )  ,  which( colnames(indata.original) %in% colnames(indata)[ samples.indata[[2]] ] )  )		

		n = lapply( samples.indata.original, length )
				
		
		indata = cbind( indata, 	rowMeans(indata.original[,samples.indata.original[[1]],drop=F]) - rowMeans(indata.original[,samples.indata.original[[2]],drop=F])  )	
		metadata = cbind( metadata, 	rowMeans(metadata[,samples.indata[[1]],drop=F]) - rowMeans(metadata[,samples.indata[[2]],drop=F])  )
	
		
		sd.shrink.1 = get.SD.estimate(data=indata.original, samples=samples.indata.original[[1]] )
		sd.shrink.2 = get.SD.estimate(data=indata.original, samples=samples.indata.original[[2]] )		
		
		
		t.g.m[,d] = rowMeans(indata.original[,samples.indata.original[[1]],drop=F]) - rowMeans(indata.original[,samples.indata.original[[2]],drop=F]) 
		t.g.m[,d] = t.g.m[,d] / sqrt( sd.shrink.1 / n[[1]] + sd.shrink.2 / n[[2]] )


		suppressWarnings({ try.res = try( { fdrtool.result = fdrtool( t.g.m[,d], verbose=F, plot=F ) }, silent=T ) })
		
		if( class(try.res) != "try-error" )
		{								
			p.g.m[,d] = fdrtool.result$pval
			fdr.g.m[,d] = fdrtool.result$lfdr
			Fdr.g.m[,d] = fdrtool.result$qval
		
			n.0.m[ d ] = fdrtool.result$param[1,"eta0"]
			perc.DE.m[ d ] = 1 - n.0.m[ d ]
			
		} else 
		{					
			p.g.m[,d] = order( indata[,d] ) / nrow(indata)
			fdr.g.m[,d] = p.g.m[,d]
			Fdr.g.m[,d] = p.g.m[,d]
			
			n.0.m[ d ] = 0.5
			perc.DE.m[ d ] = 1 - n.0.m[ d ]
		}
		
		
		delta.e.g.m = rowMeans(indata.original[,samples.indata.original[[1]],drop=F]) - rowMeans(indata.original[,samples.indata.original[[2]],drop=F])
		w.g.m = ( delta.e.g.m - min(delta.e.g.m) ) / ( max(delta.e.g.m) - min(delta.e.g.m) )
		WAD.g.m[,d] = w.g.m * delta.e.g.m		
		
	}
	
	indata = indata[ , -c(1:length(group.labels)), drop=F ]	
	colnames(indata) = names(preferences$differences.list)
	
	metadata = metadata[ , -c(1:length(group.labels)), drop=F ]	
	colnames(metadata) = names(preferences$differences.list)
	
	group.labels = names(preferences$differences.list)
	names(group.labels) = names(preferences$differences.list)
	group.colors = rep("gray20",length(preferences$differences.list))
	names(group.colors) = names(preferences$differences.list)

	
	
	
	
 	
 	output.paths = c(	"CSV" = paste( files.name, "- Results/Summary Sheets - Differences/CSV Sheets" ),
 										"Summary Sheets Samples"= paste( files.name, "- Results/Summary Sheets - Differences/Reports" ) )
 	
 	
	capture.output({	
	 	source("lib/source/Detect Spots Samples.r", local=T);
	 	
	 	if( preferences$geneset.analysis )
	 	{
			if( ncol(t.g.m) == 1 ) t.g.m = cbind(t.g.m,t.g.m)			# crack for by command, which requires >=2 columns
			source("lib/source/Geneset Statistic Samples.r", local=T);
	 	}
	 	
	 	source("lib/source/Gene Lists.r", local=T);	
	 	
	 	source("lib/source/Summary Sheets Samples.r", local=T);
	}) 	
 	 	
}
 

if( length(preferences$differences.list) > 0 )		 	
	
 	suppressWarnings( perform.differences.analyses() )
 
