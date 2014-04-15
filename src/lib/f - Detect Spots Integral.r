

	metadata.scaled = apply( metadata, 2, function(x) (x-min(x))/(max(x)-min(x)) )
	
	
	
##### Overexpression Spots ######
	
	## extract sample modules ##
	
	sample.spot.list = list()
	sample.spot.core.list = list()
	n.sample.modules = 0
	
	
	
	for( m in 1:ncol(indata) )
	{
		
		# define bigger core regions
		
		core = matrix( NA, preferences$dim.som1, preferences$dim.som1 )	
		core[ which( metadata.scaled[,m] > preferences$summary.spot.threshold ) ] = -1
		
		spot.i = 0
		while( nrow( which( core == -1, arr.ind=T ) ) > 0 )
		{
			start.pix = which( core == -1, arr.ind=T )[1,]
			
			spot.i = spot.i + 1
			core = col.pix( core, start.pix[1], start.pix[2], spot.i )
		}		
		
		
		# shrink each separate region to core size
		
		for( s.i in 1:max(core,na.rm=T) )
		{
# 			unexpressed.core.metagenes = which(core==s.i)[ which( metadata.scaled[which(core==s.i),m] < preferences$summary.spot.threshold ) ]
# 			if( length(unexpressed.core.metagenes) > 0 ) core[unexpressed.core.metagenes] = NA
						
			if( sum( core == s.i, na.rm=T ) > preferences$summary.spot.core )
			{
				core.metagenes = which( core == s.i )
				o = order( metadata[core.metagenes,m], decreasing=T )[1:preferences$summary.spot.core]
				core[ setdiff( core.metagenes, core.metagenes[o] ) ] = NA
			}			
		}
		
		core[which(!is.na(core))] = -1
		spot.i = 0
		while( nrow( which( core == -1, arr.ind=T ) ) > 0 )
		{
			start.pix = which( core == -1, arr.ind=T )[1,]
			
			spot.i = spot.i + 1
			core = col.pix( core, start.pix[1], start.pix[2], spot.i )
		}	
		
		
		
		# define spot area around cores
		
		for( s.i in 1:max(core,na.rm=T) )
		{
			n.sample.modules = n.sample.modules + 1 
			
			sample.spot.core.list[[n.sample.modules]] = matrix( NA, preferences$dim.som1, preferences$dim.som1 )
			sample.spot.core.list[[n.sample.modules]][ which( core == s.i ) ] = metadata[which( core == s.i ),m]
			
			
			
			spot = matrix( NA, preferences$dim.som1, preferences$dim.som1 )	
			spot[  which( metadata.scaled[,m] > preferences$summary.spot.threshold ) ] = -1								
			
			start.pix = which( !is.na( sample.spot.core.list[[n.sample.modules]] ) , arr.ind=T )
			start.pix = start.pix[ which.max( sample.spot.core.list[[n.sample.modules]][start.pix] ), ]
			
			spot = col.pix( spot, start.pix[1], start.pix[2], 1 )
			
			
			sample.spot.list[[n.sample.modules]] = matrix( NA, preferences$dim.som1, preferences$dim.som1 )
			sample.spot.list[[n.sample.modules]][ which( spot == 1 ) ] = metadata.scaled[which( spot == 1 ),m]
			
		}
		
		
		
	}
	
	
	
	# filter	
	
	remove = c()
	for( i in 1:n.sample.modules)	
	{
		if(	sum( !is.na( sample.spot.list[[i]] ) ) <= 1 )
		{
			# empty, i.e. does not exceed threshold -> remove
			
			remove = c( remove, i )
			
		} else
			if( sum( !is.na( sample.spot.list[[i]] ) ) < sum( !is.na( sample.spot.core.list[[i]] ) ) )
			{
				# core larger than spot -> shrink core
				
				sample.spot.core.list[[i]] = sample.spot.list[[i]]
			}
		
	}
	if( length(remove) > 0 )
	{
		sample.spot.list = sample.spot.list[-remove]
		sample.spot.core.list = sample.spot.core.list[-remove]
		n.sample.modules = length(sample.spot.list)
	}
	
	
	o = order( sapply( sample.spot.core.list,function(x) mean(x,na.rm=T) ), decreasing=T )
	sample.spot.list = sample.spot.list[o]
	sample.spot.core.list = sample.spot.core.list[o]
	
	
	
	
	
	
	
	## merge overlapping sample cores ##
	
	merged=T
	while(merged==T)
	{
		merged=F
		
		i = 1
		while( i < length(sample.spot.list) )	
		{
			
			j = i+1
			while( j <= length(sample.spot.list) )
			{
				
				core1 = which( !is.na( sample.spot.core.list[[i]] ) )
				core2 = which( !is.na( sample.spot.core.list[[j]] ) )
				
				if( any(core1 %in% core2) )
				{				
					
					# merge cores					
					if( length(setdiff(core1,core2)) > 0 )
					{		
						sample.spot.core.list[[i]][ setdiff(core2,core1) ] = 	
							sample.spot.core.list[[j]][ setdiff(core2,core1) ]
					}
					
					# merge spots
					spot1 = which( !is.na( sample.spot.list[[i]] ) )
					spot2 = which( !is.na( sample.spot.list[[j]] ) )
					
					if( length(setdiff(spot2,spot1)) > 0 )
					{
						sample.spot.list[[i]][ setdiff(spot2,spot1) ] = 	
							sample.spot.list[[j]][ setdiff(spot2,spot1) ]
					}
					
					# remove j
					sample.spot.list = sample.spot.list[-j]
					sample.spot.core.list = sample.spot.core.list[-j]
					
					merged = T
				} else
				{
					j = j + 1
				}
				
			}
			
			i = i + 1
			
		}
		
	}
	
	o = order( sapply( sample.spot.core.list,function(x) mean(x,na.rm=T) ), decreasing=T )
	sample.spot.list = sample.spot.list[o]
	sample.spot.core.list = sample.spot.core.list[o]
	
	
	
	
	
	## shrinking overlapping spots ##
	
	if( length(sample.spot.list) > 1 )	
	{ 
		for( i in 1:(length(sample.spot.list)-1) )
		{
			for( j in (i+1):length(sample.spot.list) )
			{
				
				spot1 = which( !is.na( sample.spot.list[[i]] ) )
				spot2 = which( !is.na( sample.spot.list[[j]] ) )
				
				
				if( any(spot1 %in% spot2) )
				{
					
					spot12intersect = which(spot1 %in% spot2)
					
					spot1 = which( !is.na( sample.spot.list[[i]] ), arr.ind=T )
					spot2 = which( !is.na( sample.spot.list[[j]] ), arr.ind=T )
					
					spot12intersect = spot1[ spot12intersect, ,drop=F]
					
					core1.center = colMeans( apply( which( !is.na( sample.spot.core.list[[i]] ), arr.ind=T ), 2, range ) )
					core2.center = colMeans( apply( which( !is.na( sample.spot.core.list[[j]] ), arr.ind=T ), 2, range ) )
					
					
					spot12assoc = apply( spot12intersect, 1, function(x)
					{								
						which.min(c(  sum( ( core1.center - x ) ^ 2 ), sum( ( core2.center - x ) ^ 2 )  ))
					} )
					
					
					sample.spot.list[[j]][  spot12intersect[ which(spot12assoc==1), ,drop=F]  ] = NA
					sample.spot.list[[i]][  spot12intersect[ which(spot12assoc==2), ,drop=F]  ] = NA			
					
				}
				
			}			
		}
	}
	
	
	## define overexpression spots ##
	
	
	GS.infos.overexpression = list()
	GS.infos.overexpression$overview.map = apply( apply( metadata, 2, function(x){ ( x - min(x) ) / ( max(x) - min(x) ) } ), 1, max )
	GS.infos.overexpression$overview.mask = rep( NA, preferences$dim.som1^2 )
	GS.infos.overexpression$filtered = F
	GS.infos.overexpression$spots = list()
	
	for( i in 1:length(sample.spot.list) )
	{
		GS.infos.overexpression$overview.mask[ which(!is.na(sample.spot.list[[i]])) ] = i
		
		spot.metagenes = which( !is.na(sample.spot.list[[i]]) )
		spot.genes = rownames( indata )[ which( som.nodes %in% spot.metagenes ) ]
		
		GS.infos.overexpression$spots[[ LETTERS[i] ]] = list()
		
		GS.infos.overexpression$spots[[ LETTERS[i] ]]$metagenes = spot.metagenes
		GS.infos.overexpression$spots[[ LETTERS[i] ]]$genes = spot.genes
		
		GS.infos.overexpression$spots[[ LETTERS[i] ]]$mask = rep( NA, (preferences$dim.som1*preferences$dim.som1) )
		GS.infos.overexpression$spots[[ LETTERS[i] ]]$mask[ spot.metagenes ] = 1
		
		GS.infos.overexpression$spots[[ LETTERS[i] ]]$position = colMeans( apply( som.result$code.sum[ spot.metagenes, 1:2 ]+1, 2, range ) )
		
		GS.infos.overexpression$spots[[ LETTERS[i] ]]$beta.statistic = get.beta.statistic( set.data=metadata[ GS.infos.overexpression$spots[[ LETTERS[i] ]]$metagenes,,drop=F ], weights=som.result$code.sum[ GS.infos.overexpression$spots[[ LETTERS[i] ]]$metagenes, ]$nobs )
		
	}
	
	
	o = order( sapply( GS.infos.overexpression$spots, function(x)
	{
		mean.spot.metagene = apply( metadata[ x$metagenes, ,drop=F], 2, mean )
		return(  which.max( mean.spot.metagene )  )
	} ) )	
	
	GS.infos.overexpression$spots = GS.infos.overexpression$spots[o]
	names(GS.infos.overexpression$spots) = LETTERS[ 1:length(GS.infos.overexpression$spots) ]
	
	GS.infos.overexpression$overview.mask[ !is.na(GS.infos.overexpression$overview.mask) ] = match( GS.infos.overexpression$overview.mask[ !is.na(GS.infos.overexpression$overview.mask) ], o )
	
	
	GS.infos.overexpression$spotdata = t( sapply( GS.infos.overexpression$spots, function(x) if( length( x$genes > 0 ) )	colMeans(indata[x$genes,,drop=F]) else rep( 0, ncol(indata) ) ) )
	colnames(GS.infos.overexpression$spotdata) = colnames(indata)	
	
	

	
	

##### Underexpression Spots ######
	
	## extract sample modules ##
	
	sample.spot.list = list()
	sample.spot.core.list = list()
	n.sample.modules = 0
	
	
	
	for( m in 1:ncol(indata) )
	{
		
		# define bigger core regions
		
		core = matrix( NA, preferences$dim.som1, preferences$dim.som1 )	
		core[ which( metadata.scaled[,m] < 1-preferences$summary.spot.threshold ) ] = -1
		
		spot.i = 0
		while( nrow( which( core == -1, arr.ind=T ) ) > 0 )
		{
			start.pix = which( core == -1, arr.ind=T )[1,]
			
			spot.i = spot.i + 1
			core = col.pix( core, start.pix[1], start.pix[2], spot.i )
		}		
		
		
		# shrink each separate region to core size
		
		for( s.i in 1:max(core,na.rm=T) )
		{
			if( sum( core == s.i, na.rm=T ) > preferences$summary.spot.core )
			{
				core.metagenes = which( core == s.i )
				o = order( metadata[core.metagenes,m], decreasing=F )[1:preferences$summary.spot.core]
				core[ setdiff( core.metagenes, core.metagenes[o] ) ] = NA
			}			
		}
		
		core[which(!is.na(core))] = -1
		spot.i = 0
		while( nrow( which( core == -1, arr.ind=T ) ) > 0 )
		{
			start.pix = which( core == -1, arr.ind=T )[1,]
			
			spot.i = spot.i + 1
			core = col.pix( core, start.pix[1], start.pix[2], spot.i )
		}	
		
		
		
		# define spot area around cores
		
		for( s.i in 1:max(core,na.rm=T) )
		{
			n.sample.modules = n.sample.modules + 1 
			
			sample.spot.core.list[[n.sample.modules]] = matrix( NA, preferences$dim.som1, preferences$dim.som1 )
			sample.spot.core.list[[n.sample.modules]][ which( core == s.i ) ] = metadata[which( core == s.i ),m]
			
			
			
			spot = matrix( NA, preferences$dim.som1, preferences$dim.som1 )	
			spot[  which( metadata.scaled[,m] < 1-preferences$summary.spot.threshold ) ] = -1								
			
			start.pix = which( !is.na( sample.spot.core.list[[n.sample.modules]] ) , arr.ind=T )
			start.pix = start.pix[ which.max( sample.spot.core.list[[n.sample.modules]][start.pix] ), ]
			
			spot = col.pix( spot, start.pix[1], start.pix[2], 1 )
			
			
			sample.spot.list[[n.sample.modules]] = matrix( NA, preferences$dim.som1, preferences$dim.som1 )
			sample.spot.list[[n.sample.modules]][ which( spot == 1 ) ] = metadata.scaled[which( spot == 1 ),m]
			
		}
		
		
		
	}
	
	
	
	# filter	
	
	remove = c()
	for( i in 1:n.sample.modules)	
	{
		if(	sum( !is.na( sample.spot.list[[i]] ) ) <= 1 )
		{
			# empty, i.e. does not exceed threshold -> remove
			
			remove = c( remove, i )
			
		} else
			if( sum( !is.na( sample.spot.list[[i]] ) ) < sum( !is.na( sample.spot.core.list[[i]] ) ) )
			{
				# core larger than spot -> shrink core
				
				sample.spot.core.list[[i]] = sample.spot.list[[i]]
			}
		
	}
	if( length(remove) > 0 )
	{
		sample.spot.list = sample.spot.list[-remove]
		sample.spot.core.list = sample.spot.core.list[-remove]
		n.sample.modules = length(sample.spot.list)
	}
	
	
	o = order( sapply( sample.spot.core.list,function(x) mean(x,na.rm=T) ), decreasing=F )
	sample.spot.list = sample.spot.list[o]
	sample.spot.core.list = sample.spot.core.list[o]
	
	
	
	
	
	
	
	## merge overlapping sample cores ##
	
	merged=T
	while(merged==T)
	{
		merged=F
		
		i = 1
		while( i < length(sample.spot.list) )	
		{
			
			j = i+1
			while( j <= length(sample.spot.list) )
			{
				
				core1 = which( !is.na( sample.spot.core.list[[i]] ) )
				core2 = which( !is.na( sample.spot.core.list[[j]] ) )
				
				if( any(core1 %in% core2) )
				{				
					
					# merge cores					
					if( length(setdiff(core1,core2)) > 0 )
					{		
						sample.spot.core.list[[i]][ setdiff(core2,core1) ] = 	
							sample.spot.core.list[[j]][ setdiff(core2,core1) ]
					}
					
					# merge spots
					spot1 = which( !is.na( sample.spot.list[[i]] ) )
					spot2 = which( !is.na( sample.spot.list[[j]] ) )
					
					if( length(setdiff(spot2,spot1)) > 0 )
					{
						sample.spot.list[[i]][ setdiff(spot2,spot1) ] = 	
							sample.spot.list[[j]][ setdiff(spot2,spot1) ]
					}
					
					# remove j
					sample.spot.list = sample.spot.list[-j]
					sample.spot.core.list = sample.spot.core.list[-j]
					
					merged = T
				} else
				{
					j = j + 1
				}
				
			}
			
			i = i + 1
			
		}
		
	}
	
	o = order( sapply( sample.spot.core.list,function(x) mean(x,na.rm=T) ), decreasing=F )
	sample.spot.list = sample.spot.list[o]
	sample.spot.core.list = sample.spot.core.list[o]
	
	
	
	
	
	## shrinking overlapping spots ##
	
	if( length(sample.spot.list) > 1 )	
	{ 
		for( i in 1:(length(sample.spot.list)-1) )
		{
			for( j in (i+1):length(sample.spot.list) )
			{
				
				spot1 = which( !is.na( sample.spot.list[[i]] ) )
				spot2 = which( !is.na( sample.spot.list[[j]] ) )
				
				
				if( any(spot1 %in% spot2) )
				{
					
					spot12intersect = which(spot1 %in% spot2)
					
					spot1 = which( !is.na( sample.spot.list[[i]] ), arr.ind=T )
					spot2 = which( !is.na( sample.spot.list[[j]] ), arr.ind=T )
					
					spot12intersect = spot1[ spot12intersect, ,drop=F]
					
					core1.center = colMeans( apply( which( !is.na( sample.spot.core.list[[i]] ), arr.ind=T ), 2, range ) )
					core2.center = colMeans( apply( which( !is.na( sample.spot.core.list[[j]] ), arr.ind=T ), 2, range ) )
					
					
					spot12assoc = apply( spot12intersect, 1, function(x)
					{								
						which.min(c(  sum( ( core1.center - x ) ^ 2 ), sum( ( core2.center - x ) ^ 2 )  ))
					} )
					
					
					sample.spot.list[[j]][  spot12intersect[ which(spot12assoc==1), ,drop=F]  ] = NA
					sample.spot.list[[i]][  spot12intersect[ which(spot12assoc==2), ,drop=F]  ] = NA			
					
				}
				
			}
		}
	}
	
	
	## define underexpression spots ##
	
	
	GS.infos.underexpression = list()
	GS.infos.underexpression$overview.map = apply( apply( metadata, 2, function(x){ ( x - min(x) ) / ( max(x) - min(x) ) } ), 1, min )
	GS.infos.underexpression$overview.mask = rep( NA, preferences$dim.som1^2 )
	GS.infos.underexpression$filtered = F
	GS.infos.underexpression$spots = list()
	
	for( i in 1:length(sample.spot.list) )
	{
		GS.infos.underexpression$overview.mask[ which(!is.na(sample.spot.list[[i]])) ] = i
		
		spot.metagenes = which( !is.na(sample.spot.list[[i]]) )
		spot.genes = rownames( indata )[ which( som.nodes %in% spot.metagenes ) ]
		
		GS.infos.underexpression$spots[[ letters[i] ]] = list()
		
		GS.infos.underexpression$spots[[ letters[i] ]]$metagenes = spot.metagenes
		GS.infos.underexpression$spots[[ letters[i] ]]$genes = spot.genes
		
		GS.infos.underexpression$spots[[ letters[i] ]]$mask = rep( NA, (preferences$dim.som1*preferences$dim.som1) )
		GS.infos.underexpression$spots[[ letters[i] ]]$mask[ spot.metagenes ] = 1
		
		GS.infos.underexpression$spots[[ letters[i] ]]$position = colMeans( apply( som.result$code.sum[ spot.metagenes, 1:2 ]+1, 2, range ) )
		
		GS.infos.underexpression$spots[[ letters[i] ]]$beta.statistic = get.beta.statistic( set.data=metadata[ GS.infos.underexpression$spots[[ letters[i] ]]$metagenes,,drop=F ], weights=som.result$code.sum[ GS.infos.underexpression$spots[[ letters[i] ]]$metagenes, ]$nobs )
		
	}
	
	
	o = order( sapply( GS.infos.underexpression$spots, function(x)
	{
		mean.spot.metagene = apply( metadata[ x$metagenes, ,drop=F], 2, mean )
		return(  which.min( mean.spot.metagene )  )
	} ) )	
	
	GS.infos.underexpression$spots = GS.infos.underexpression$spots[o]
	names(GS.infos.underexpression$spots) = letters[ 1:length(GS.infos.underexpression$spots) ]
	
	GS.infos.underexpression$overview.mask[ !is.na(GS.infos.underexpression$overview.mask) ] = match( GS.infos.underexpression$overview.mask[ !is.na(GS.infos.underexpression$overview.mask) ], o )
	
	
	GS.infos.underexpression$spotdata = t( sapply( GS.infos.underexpression$spots, function(x) if( length( x$genes > 0 ) )	colMeans(indata[x$genes,,drop=F]) else rep( 0, ncol(indata) ) ) )
	colnames(GS.infos.underexpression$spotdata) = colnames(indata)	
	
	
	 

	
		
	
# 	
# 	
# 
# 	##### Positive Peaks ######
# 
# 	
# 	peaks = apply( metadata, 1, max )
# 
# 
# 	GS.infos.positivepeaks = list()
# 	GS.infos.positivepeaks$overview.map = peaks
# 	GS.infos.positivepeaks$overview.mask = rep( NA, preferences$dim.som1^2 )
# 	GS.infos.positivepeaks$filtered = F
# 	GS.infos.positivepeaks$spots = list()	
# 
# 
# 	peaks = matrix( peaks, preferences$dim.som1, preferences$dim.som1 )
# 	peaks[ which( peaks < quantile( peaks, preferences$summary.spot.threshold ) ) ] = NA
# 
# 	e.cluster = matrix( NA, preferences$dim.som1, preferences$dim.som1 )
# 	e.cluster[ which( !is.na(peaks) ) ] = -1
# 	count.cluster = 1	
# 
# 
# 	while( any( !is.na( peaks ) ) )
# 	{
# 		
# 		start.pix = which( peaks == max(peaks,na.rm=T), arr.ind=T )[1,]
# 
# 		e.cluster = col.pix( e.cluster, start.pix[1], start.pix[2], count.cluster )
# 
# 		peaks[ which( e.cluster == count.cluster ) ] = NA
# 
# 		count.cluster = count.cluster + 1
# 	}
# 
# 
# 
# 	e.cluster = as.vector( e.cluster )
# 
# 	for( i in 1:length( na.omit( unique( e.cluster ) ) ) )
# 	{
# 		nodes = which( e.cluster == i )
# 
# 
# 		geneset.genes = rownames( indata )[ which( som.nodes %in% nodes ) ]
# 
# 		
# 		GS.infos.positivepeaks$overview.mask[ as.numeric( nodes ) ] = i
# 		
# 		GS.infos.positivepeaks$spots[[ LETTERS[i] ]] = list()
# 
# 		GS.infos.positivepeaks$spots[[ LETTERS[i] ]]$metagenes = as.numeric( nodes )
# 		GS.infos.positivepeaks$spots[[ LETTERS[i] ]]$genes = geneset.genes
# 
# 		GS.infos.positivepeaks$spots[[ LETTERS[i] ]]$mask = rep( NA, (preferences$dim.som1*preferences$dim.som1) )
# 		GS.infos.positivepeaks$spots[[ LETTERS[i] ]]$mask[ as.numeric( nodes ) ] = 1
# 
# 		GS.infos.positivepeaks$spots[[ LETTERS[i] ]]$position = apply( apply( som.result$code.sum[ nodes, 1:2 ], 2, range ), 2, mean ) + 0.5
# 		
# 		GS.infos.positivepeaks$spots[[ LETTERS[i] ]]$beta.statistic = get.beta.statistic( set.data=metadata[ GS.infos.positivepeaks$spots[[ LETTERS[i] ]]$metagenes,,drop=F ], weights=som.result$code.sum[ GS.infos.positivepeaks$spots[[ LETTERS[i] ]]$metagenes, ]$nobs )
# 	}
# 
# 	GS.infos.positivepeaks$spotdata = t( sapply( GS.infos.positivepeaks$spots, function(x) if( length( x$genes > 0 ) )	colMeans(indata[x$genes,,drop=F]) else rep( 0, ncol(indata) ) ) )
# 	colnames(GS.infos.positivepeaks$spotdata) = colnames(indata)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 	##### Negative Peaks ######
# 
# 	
# 	peaks = apply( metadata, 1, min )
# 
# 
# 	GS.infos.negativepeaks = list()
# 	GS.infos.negativepeaks$overview.map = peaks
# 	GS.infos.negativepeaks$overview.mask = rep( NA, preferences$dim.som1^2 )
# 	GS.infos.negativepeaks$filtered = F
# 	GS.infos.negativepeaks$spots = list()
# 
# 
# 
# 
# 	peaks = matrix( peaks, preferences$dim.som1, preferences$dim.som1 )
# 	peaks[ which( peaks > quantile( peaks, 1-preferences$summary.spot.threshold ) ) ] = NA
# 
# 	e.cluster = matrix( NA, preferences$dim.som1, preferences$dim.som1 )
# 	e.cluster[ which( !is.na(peaks) ) ] = -1
# 	count.cluster = 1	
# 
# 
# 	while( any( !is.na( peaks ) ) )
# 	{
# 		
# 		start.pix = which( peaks == min(peaks,na.rm=T), arr.ind=T )[1,]
# 
# 		e.cluster = col.pix( e.cluster, start.pix[1], start.pix[2], count.cluster )
# 
# 		peaks[ which( e.cluster == count.cluster ) ] = NA
# 
# 		count.cluster = count.cluster + 1
# 	}
# 
# 
# 
# 	e.cluster = as.vector( e.cluster )
# 
# 	for( i in 1:length( na.omit( unique( e.cluster ) ) ) )
# 	{
# 		nodes = which( e.cluster == i )
# 
# 
# 		geneset.genes = rownames( indata )[ which( som.nodes %in% nodes ) ]
# 
# 		
# 		GS.infos.negativepeaks$overview.mask[ as.numeric( nodes ) ] = i		
# 
# 		GS.infos.negativepeaks$spots[[ letters[i] ]] = list()
# 
# 		GS.infos.negativepeaks$spots[[ letters[i] ]]$metagenes = as.numeric( nodes )
# 		GS.infos.negativepeaks$spots[[ letters[i] ]]$genes = geneset.genes
# 
# 		GS.infos.negativepeaks$spots[[ letters[i] ]]$mask = rep( NA, (preferences$dim.som1*preferences$dim.som1) )
# 		GS.infos.negativepeaks$spots[[ letters[i] ]]$mask[ as.numeric( nodes ) ] = 1
# 
# 		GS.infos.negativepeaks$spots[[ letters[i] ]]$position = apply( apply( som.result$code.sum[ nodes, 1:2 ], 2, range ), 2, mean ) + 0.5		
# 		
# 		GS.infos.negativepeaks$spots[[ letters[i] ]]$beta.statistic = get.beta.statistic( set.data=metadata[ GS.infos.negativepeaks$spots[[ letters[i] ]]$metagenes,,drop=F ], weights=som.result$code.sum[ GS.infos.negativepeaks$spots[[ letters[i] ]]$metagenes, ]$nobs )
# 	}
# 
# 	GS.infos.negativepeaks$spotdata = t( sapply( GS.infos.negativepeaks$spots, function(x) if( length( x$genes > 0 ) )	colMeans(indata[x$genes,,drop=F]) else rep( 0, ncol(indata) ) ) )
# 	colnames(GS.infos.negativepeaks$spotdata) = colnames(indata)
# 
# 













	##### Correlation Cluster ######


	GS.infos.correlation = list()
	GS.infos.correlation$overview.map = NA	
	GS.infos.correlation$overview.mask = rep( NA, preferences$dim.som1^2 )
	GS.infos.correlation$filtered = F
	GS.infos.correlation$spots = list()



	c.map = cor( t( metadata ) )
	diag( c.map ) = NA
	rownames( c.map ) = c( 1:( preferences$dim.som1*preferences$dim.som1 ) )
	colnames( c.map ) = c( 1:( preferences$dim.som1*preferences$dim.som1 ) )

	c.cluster = rep( NA, (preferences$dim.som1*preferences$dim.som1) )
	names( c.cluster ) = c( 1:( preferences$dim.som1*preferences$dim.som1 ) )



	count.cluster = 1


	while( is.matrix(c.map) && nrow(c.map) > 0 && ncol(c.map) > 0 && count.cluster <= 30 )
	{
		
		
		start.node = rownames(c.map)[ which( c.map == max( c.map, na.rm=T ), arr.ind=T )[1,1] ]

		cluster = names( which( c.map[ start.node, ] > 0.90 ) )
		cluster = c( start.node, cluster )


		if( length(cluster) >= preferences$dim.som1 / 2 )
		{
			c.cluster[ cluster ] = count.cluster

			geneset.genes = rownames( indata )[ which( som.nodes %in% as.numeric( cluster ) ) ]

		
			GS.infos.correlation$overview.mask[ as.numeric( cluster ) ] = count.cluster
			
			GS.infos.correlation$spots[[ count.cluster ]] = list()

			GS.infos.correlation$spots[[ count.cluster ]]$metagenes = as.numeric( cluster )
			GS.infos.correlation$spots[[ count.cluster ]]$genes = geneset.genes

			GS.infos.correlation$spots[[ count.cluster ]]$mask = rep( NA, (preferences$dim.som1*preferences$dim.som1) )
			GS.infos.correlation$spots[[ count.cluster ]]$mask[ as.numeric( cluster) ] = 1

			GS.infos.correlation$spots[[ count.cluster ]]$position = apply( apply( som.result$code.sum[ cluster, 1:2 ], 2, range ), 2, mean ) + 0.5

			GS.infos.correlation$spots[[ count.cluster ]]$beta.statistic = get.beta.statistic( set.data=metadata[ GS.infos.correlation$spots[[ count.cluster ]]$metagenes,,drop=F ], weights=som.result$code.sum[ GS.infos.correlation$spots[[ count.cluster ]]$metagenes, ]$nobs )
			
			count.cluster = count.cluster + 1
		}

		c.map = c.map[ -which( rownames(c.map)%in%cluster), -which( colnames(c.map)%in%cluster) ]


	}
	GS.infos.correlation$overview.map = c.cluster


	o = order( sapply( GS.infos.correlation$spots, function(x)
	{
		mean.spot.metagene = apply( metadata[ x$metagenes, ,drop=F], 2, mean )
		return(  which.max( mean.spot.metagene )  )
	} ) )

	GS.infos.correlation$spots = GS.infos.correlation$spots[o]	
	names(GS.infos.correlation$spots) = LETTERS[ 1:length(GS.infos.correlation$spots) ]

	GS.infos.correlation$overview.mask[ !is.na(GS.infos.correlation$overview.mask) ] = match( GS.infos.correlation$overview.mask[ !is.na(GS.infos.correlation$overview.mask) ], o )


	GS.infos.correlation$spotdata = t( sapply( GS.infos.correlation$spots, function(x) if( length( x$genes > 0 ) )	colMeans(indata[x$genes,,drop=F]) else rep( 0, ncol(indata) ) ) )
	colnames(GS.infos.correlation$spotdata) = colnames(indata)	
	




	##### K-Means Clustering #####

	
	n.cluster = preferences$dim.som1 / 2
	
	prototypes = metadata[ round( seq( 1, preferences$dim.som1^2, length.out=n.cluster ) ), ]
		
	res = kmeans( metadata, prototypes )
	


	GS.infos.kmeans = list()
	GS.infos.kmeans$overview.map = matrix( res$cluster, preferences$dim.som1, preferences$dim.som1 )
	GS.infos.kmeans$overview.mask = rep( NA, preferences$dim.som1^2 )
	GS.infos.kmeans$filtered = F
	GS.infos.kmeans$spots = list()
	

	for( i in 1:n.cluster )
	{
		nodes = which( res$cluster == i )


		geneset.genes = rownames( indata )[ which( som.nodes %in% nodes ) ]

		
		GS.infos.kmeans$overview.mask[ as.numeric( nodes ) ] = i		

		GS.infos.kmeans$spots[[ LETTERS[i] ]] = list()

		GS.infos.kmeans$spots[[ LETTERS[i] ]]$metagenes = as.numeric( nodes )
		GS.infos.kmeans$spots[[ LETTERS[i] ]]$genes = geneset.genes

		GS.infos.kmeans$spots[[ LETTERS[i] ]]$mask = rep( NA, (preferences$dim.som1*preferences$dim.som1) )
		GS.infos.kmeans$spots[[ LETTERS[i] ]]$mask[ as.numeric( nodes ) ] = 1

		GS.infos.kmeans$spots[[ LETTERS[i] ]]$position = apply( apply( som.result$code.sum[ nodes, 1:2 ], 2, range ), 2, mean ) + 0.5
		
		GS.infos.kmeans$spots[[ LETTERS[i] ]]$beta.statistic = get.beta.statistic( set.data=metadata[ GS.infos.kmeans$spots[[ LETTERS[i] ]]$metagenes,,drop=F ], weights=som.result$code.sum[ GS.infos.kmeans$spots[[ LETTERS[i] ]]$metagenes, ]$nobs )
	}


	o = order( sapply( GS.infos.kmeans$spots, function(x)
	{
		mean.spot.metagene = apply( metadata[ x$metagenes, ,drop=F], 2, mean )
		return(  which.max( mean.spot.metagene )  )
	} ) )

	GS.infos.kmeans$spots = GS.infos.kmeans$spots[o]	
	names(GS.infos.kmeans$spots) = LETTERS[ 1:length(GS.infos.kmeans$spots) ]
	
	GS.infos.kmeans$overview.mask[ !is.na(GS.infos.kmeans$overview.mask) ] = match( GS.infos.kmeans$overview.mask[ !is.na(GS.infos.kmeans$overview.mask) ], o )

	
	GS.infos.kmeans$spotdata = t( sapply( GS.infos.kmeans$spots, function(x) if( length( x$genes > 0 ) )	colMeans(indata[x$genes,,drop=F]) else rep( 0, ncol(indata) ) ) )
	colnames(GS.infos.kmeans$spotdata) = colnames(indata)		
	
	
	
	
	
	
	##### Group Spots ######
	
	
	if( length(unique( group.labels ) ) > 1 )	
	{
		
		group.metadata.scaled = apply( group.metadata, 2, function(x) (x-min(x))/(max(x)-min(x)) )

		
			
		## extract sample modules ##
		
		sample.spot.list = list()
		sample.spot.core.list = list()
		n.sample.modules = 0
		
		
		
		for( m in 1:ncol(group.metadata) )
		{
			
			# define bigger core regions
			
			core = matrix( NA, preferences$dim.som1, preferences$dim.som1 )	
			core[ which( group.metadata.scaled[,m] > preferences$group.spot.threshold ) ] = -1
			
			spot.i = 0
			while( nrow( which( core == -1, arr.ind=T ) ) > 0 )
			{
				start.pix = which( core == -1, arr.ind=T )[1,]
				
				spot.i = spot.i + 1
				core = col.pix( core, start.pix[1], start.pix[2], spot.i )
			}		
			
			
			# shrink each separate region to core size
			
			for( s.i in 1:max(core,na.rm=T) )
			{
				# 			unexpressed.core.metagenes = which(core==s.i)[ which( group.metadata.scaled[which(core==s.i),m] < preferences$group.spot.threshold ) ]
				# 			if( length(unexpressed.core.metagenes) > 0 ) core[unexpressed.core.metagenes] = NA
				
				if( sum( core == s.i, na.rm=T ) > preferences$group.spot.core )
				{
					core.metagenes = which( core == s.i )
					o = order( group.metadata[core.metagenes,m], decreasing=T )[1:preferences$group.spot.core]
					core[ setdiff( core.metagenes, core.metagenes[o] ) ] = NA
				}			
			}
			
			core[which(!is.na(core))] = -1
			spot.i = 0
			while( nrow( which( core == -1, arr.ind=T ) ) > 0 )
			{
				start.pix = which( core == -1, arr.ind=T )[1,]
				
				spot.i = spot.i + 1
				core = col.pix( core, start.pix[1], start.pix[2], spot.i )
			}	
			
			
			
			# define spot area around cores
			
			for( s.i in 1:max(core,na.rm=T) )
			{
				n.sample.modules = n.sample.modules + 1 
				
				sample.spot.core.list[[n.sample.modules]] = matrix( NA, preferences$dim.som1, preferences$dim.som1 )
				sample.spot.core.list[[n.sample.modules]][ which( core == s.i ) ] = metadata[which( core == s.i ),m]
				
				
				
				spot = matrix( NA, preferences$dim.som1, preferences$dim.som1 )	
				spot[  which( group.metadata.scaled[,m] > preferences$group.spot.threshold ) ] = -1								
				
				start.pix = which( !is.na( sample.spot.core.list[[n.sample.modules]] ) , arr.ind=T )
				start.pix = start.pix[ which.max( sample.spot.core.list[[n.sample.modules]][start.pix] ), ]
				
				spot = col.pix( spot, start.pix[1], start.pix[2], 1 )
				
				
				sample.spot.list[[n.sample.modules]] = matrix( NA, preferences$dim.som1, preferences$dim.som1 )
				sample.spot.list[[n.sample.modules]][ which( spot == 1 ) ] = group.metadata.scaled[which( spot == 1 ),m]
				
			}
			
			
			
		}
		
		
		
		# filter	
		
		remove = c()
		for( i in 1:n.sample.modules)	
		{
			if(	sum( !is.na( sample.spot.list[[i]] ) ) <= 1 )
			{
				# empty, i.e. does not exceed threshold -> remove
				
				remove = c( remove, i )
				
			} else
				if( sum( !is.na( sample.spot.list[[i]] ) ) < sum( !is.na( sample.spot.core.list[[i]] ) ) )
				{
					# core larger than spot -> shrink core
					
					sample.spot.core.list[[i]] = sample.spot.list[[i]]
				}
			
		}
		if( length(remove) > 0 )
		{
			sample.spot.list = sample.spot.list[-remove]
			sample.spot.core.list = sample.spot.core.list[-remove]
			n.sample.modules = length(sample.spot.list)
		}
		
		
		o = order( sapply( sample.spot.core.list,function(x) mean(x,na.rm=T) ), decreasing=T )
		sample.spot.list = sample.spot.list[o]
		sample.spot.core.list = sample.spot.core.list[o]
		
		
		
		
		
		
		
		## merge overlapping sample cores ##
		
		merged=T
		while(merged==T)
		{
			merged=F
			
			i = 1
			while( i < length(sample.spot.list) )	
			{
				
				j = i+1
				while( j <= length(sample.spot.list) )
				{
					
					core1 = which( !is.na( sample.spot.core.list[[i]] ) )
					core2 = which( !is.na( sample.spot.core.list[[j]] ) )
					
					if( any(core1 %in% core2) )
					{				
						
						# merge cores					
						if( length(setdiff(core1,core2)) > 0 )
						{		
							sample.spot.core.list[[i]][ setdiff(core2,core1) ] = 	
								sample.spot.core.list[[j]][ setdiff(core2,core1) ]
						}
						
						# merge spots
						spot1 = which( !is.na( sample.spot.list[[i]] ) )
						spot2 = which( !is.na( sample.spot.list[[j]] ) )
						
						if( length(setdiff(spot2,spot1)) > 0 )
						{
							sample.spot.list[[i]][ setdiff(spot2,spot1) ] = 	
								sample.spot.list[[j]][ setdiff(spot2,spot1) ]
						}
						
						# remove j
						sample.spot.list = sample.spot.list[-j]
						sample.spot.core.list = sample.spot.core.list[-j]
						
						merged = T
					} else
					{
						j = j + 1
					}
					
				}
				
				i = i + 1
				
			}
			
		}
		
		o = order( sapply( sample.spot.core.list,function(x) mean(x,na.rm=T) ), decreasing=T )
		sample.spot.list = sample.spot.list[o]
		sample.spot.core.list = sample.spot.core.list[o]
		
		
		
		
		
		## shrinking overlapping spots ##
		
		if( length(sample.spot.list) > 1 )	
		{ 
			for( i in 1:(length(sample.spot.list)-1) )
			{
				for( j in (i+1):length(sample.spot.list) )
				{
					
					spot1 = which( !is.na( sample.spot.list[[i]] ) )
					spot2 = which( !is.na( sample.spot.list[[j]] ) )
					
					
					if( any(spot1 %in% spot2) )
					{
						
						spot12intersect = which(spot1 %in% spot2)
						
						spot1 = which( !is.na( sample.spot.list[[i]] ), arr.ind=T )
						spot2 = which( !is.na( sample.spot.list[[j]] ), arr.ind=T )
						
						spot12intersect = spot1[ spot12intersect, ,drop=F]
						
						core1.center = colMeans( apply( which( !is.na( sample.spot.core.list[[i]] ), arr.ind=T ), 2, range ) )
						core2.center = colMeans( apply( which( !is.na( sample.spot.core.list[[j]] ), arr.ind=T ), 2, range ) )
						
						
						spot12assoc = apply( spot12intersect, 1, function(x)
						{								
							which.min(c(  sum( ( core1.center - x ) ^ 2 ), sum( ( core2.center - x ) ^ 2 )  ))
						} )
						
						
						sample.spot.list[[j]][  spot12intersect[ which(spot12assoc==1), ,drop=F]  ] = NA
						sample.spot.list[[i]][  spot12intersect[ which(spot12assoc==2), ,drop=F]  ] = NA			
						
					}
					
				}
			}
		}
		
		
		## define overexpression spots ##
		
		
		GS.infos.group.overexpression = list()
		GS.infos.group.overexpression$overview.map = apply( apply( group.metadata, 2, function(x){ ( x - min(x) ) / ( max(x) - min(x) ) } ), 1, max )
		GS.infos.group.overexpression$overview.mask = rep( NA, preferences$dim.som1^2 )
		GS.infos.group.overexpression$filtered = F
		GS.infos.group.overexpression$spots = list()
		
		for( i in 1:length(sample.spot.list) )
		{
			GS.infos.group.overexpression$overview.mask[ which(!is.na(sample.spot.list[[i]])) ] = i
			
			spot.metagenes = which( !is.na(sample.spot.list[[i]]) )
			spot.genes = rownames( indata )[ which( som.nodes %in% spot.metagenes ) ]
			
			GS.infos.group.overexpression$spots[[ LETTERS[i] ]] = list()
			
			GS.infos.group.overexpression$spots[[ LETTERS[i] ]]$metagenes = spot.metagenes
			GS.infos.group.overexpression$spots[[ LETTERS[i] ]]$genes = spot.genes
			
			GS.infos.group.overexpression$spots[[ LETTERS[i] ]]$mask = rep( NA, (preferences$dim.som1*preferences$dim.som1) )
			GS.infos.group.overexpression$spots[[ LETTERS[i] ]]$mask[ spot.metagenes ] = 1
			
			GS.infos.group.overexpression$spots[[ LETTERS[i] ]]$position = colMeans( apply( som.result$code.sum[ spot.metagenes, 1:2 ]+1, 2, range ) )
			
			GS.infos.group.overexpression$spots[[ LETTERS[i] ]]$beta.statistic = get.beta.statistic( set.data=metadata[ GS.infos.group.overexpression$spots[[ LETTERS[i] ]]$metagenes,,drop=F ], weights=som.result$code.sum[ GS.infos.group.overexpression$spots[[ LETTERS[i] ]]$metagenes, ]$nobs )
			
		}
	
		
		
		start.spot = which.min( sapply( GS.infos.group.overexpression$spots, function(x)
		{
			mean.spot.metagene = apply( group.metadata[ x$metagenes, ,drop=F], 2, mean )	
			return(  which.max( mean.spot.metagene )  )
		} ) )		
		
		spot.arcs = sapply( GS.infos.group.overexpression$spots, function(x) 
		{
			-atan2( x$position['y']-preferences$dim.som1/2, x$position['x']-preferences$dim.som1/2 )
		})
		
		spot.arcs = spot.arcs - spot.arcs[start.spot]
		if( any(spot.arcs<0) ) spot.arcs[which(spot.arcs<0)] = spot.arcs[which(spot.arcs<0)] + (2*pi)
		
		o = order( spot.arcs )
		
		
		
		GS.infos.group.overexpression$spots = GS.infos.group.overexpression$spots[o]
		names(GS.infos.group.overexpression$spots) = LETTERS[ 1:length(GS.infos.group.overexpression$spots) ]
		
		GS.infos.group.overexpression$overview.mask[ !is.na(GS.infos.group.overexpression$overview.mask) ] = match( GS.infos.group.overexpression$overview.mask[ !is.na(GS.infos.group.overexpression$overview.mask) ], o )
		
		
		GS.infos.group.overexpression$spotdata = t( sapply( GS.infos.group.overexpression$spots, function(x) if( length( x$genes > 0 ) )	colMeans(indata[x$genes,,drop=F]) else rep( 0, ncol(indata) ) ) )
		colnames(GS.infos.group.overexpression$spotdata) = colnames(indata)	
		
		
	}	
		
		
		
		