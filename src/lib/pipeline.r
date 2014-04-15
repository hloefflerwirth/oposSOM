

	source("lib/source/Import.r")
	
	
	
	# standard preferences
	

	if( !exists("preferences") )
	{
		preferences = list()
	} 

	if( ! ( "dataset.name" %in% names(preferences) ) )		preferences$dataset.name = "Unnamed"
	
	if( ! ( "error.model" %in% names(preferences) ) )		preferences$error.model = "all.samples"
	if( ! preferences$error.model %in% c( "replicates", "all.samples", "groups" ) ) 
	{
		preferences$error.model = "all.samples"
		cat( "\n!!Unknown value of \"error.model\". Using \"all.samples\"!!\n" ); flush.console()
	}


	if( ! ( "dim.som1" %in% names(preferences) ) )			preferences$dim.som1 = 20
	if( ! ( "dim.som2" %in% names(preferences) ) )			preferences$dim.som2 = 20
	if( ! ( "training.extension" %in% names(preferences) ) )	preferences$training.extension = 1
	
	if( ! ( "rotate.som1" %in% names(preferences) ) )		preferences$rotate.som1 = 0
	
	if( ! ( "flip.som1" %in% names(preferences) ) )			preferences$flip.som1 = F

	if( ! ( "ensembl.dataset" %in% names(preferences) ) )		preferences$ensembl.dataset = ""	
	if( ! ( "ensembl.rowname.ids" %in% names(preferences) ) )	preferences$ensembl.rowname.ids = ""
	
	if( ! ( "geneset.analysis" %in% names(preferences) ) )   	preferences$geneset.analysis = F
	if( ! ( "geneset.analysis.exact" %in% names(preferences) ) )preferences$geneset.analysis.exact = F
	if( ! ( "geneset.custom.list" %in% names(preferences) ) )   preferences$geneset.custom.list= list()
	
	if( ! ( "max.parallel.cores" %in% names(preferences) ) )	preferences$max.parallel.cores = detectCores() / 2	
	
	if( ! ( "sample.spot.cutoff" %in% names(preferences) ) )	preferences$sample.spot.cutoff = 0.65
	if( ! ( "summary.spot.core" %in% names(preferences) ) )	preferences$summary.spot.core = 3
	if( ! ( "summary.spot.threshold" %in% names(preferences) ) )	preferences$summary.spot.threshold = 0.95
	if( ! ( "group.spot.core" %in% names(preferences) ) )	preferences$group.spot.core = 5
	if( ! ( "group.spot.threshold" %in% names(preferences) ) )	preferences$group.spot.threshold = 0.75
	
	if( ! ( "feature.mean.normalization" %in% names(preferences) ) )	preferences$feature.mean.normalization = T
	if( ! ( "sample.quantile.normalization" %in% names(preferences) ) )	preferences$sample.quantile.normalization = T
	
	if( ! ( "differences.list" %in% names(preferences) ) )	preferences$differences.list = list()
		
	preferences$system.info = Sys.info()
	preferences$started = format(Sys.time(), "%a %b %d %X\n" )
	
	
	
	cat( "\n\n\nStarted:", format(Sys.time(), "%a %b %d %X\n" ) )
	cat( "Setting:", preferences$dataset.name, "\n" ) 
	cat( "1SOM Dim:", preferences$dim.som1,"\n2SOM Dim:",preferences$dim.som2,"\n\n" )
	
	flush.console()
	




	
	# check input parameters/data


	if( class( indata ) != "matrix" || mode( indata ) != "numeric" || storage.mode(indata) != "numeric" )
	{
		rn = rownames( indata )
		indata = apply( indata, 2, function(x){ as.numeric( as.vector( x ) ) } )
		rownames( indata ) = rn
		storage.mode(indata) = "numeric"
		#cat( "\n!!Converted indata to numerical matrix!!\n" ); flush.console()
	}

	if( length( rownames( indata ) ) == 0 )
	{
		rownames( indata ) = as.character( 1:nrow(indata) )
		preferences$geneset.analysis = F
		cat( "\n!!No rownames found. Set them to 1,2,3,4...!!\n" ); flush.console()
	}
	
	if( length( colnames( indata ) ) == 0 )
	{
		colnames( indata ) = paste( "Sample", c( 1:ncol(indata) ) )
		cat( "\n!!No colnames found. Set them to 1,2,3,4...!!\n" ); flush.console()
	}

	if( any(duplicated(rownames(indata))) )
	{
		indata = do.call( rbind, by( indata, rownames(indata), colMeans ) )[ unique(rownames(indata)), ]
		cat( "\n!!Duplicate rownames. Averaged multiple features!!\n" ); flush.console()
	}
	
	na.rows = which( apply( apply( indata, 1, is.na ), 2, sum ) > 0 )
	if( length( na.rows ) > 0 )
	{
		indata = indata[ -na.rows, ]
		if( exists("indata.original") )
			indata.original = indata.original[ -na.rows, ]		
		
		cat( "\n!!Removed NAs from data set!!\n" ); flush.console()
	}

	
	
	
	## set up global variables
		
	
	files.name = preferences$dataset.name
	
	
	while( file.exists( paste( files.name, ".RData", sep="" ) ) )
		files.name = paste( files.name, "+", sep="" )
	
	
	
	output.paths = c( "LPE" = paste( files.name, "- Results/LPE" ),
										"CSV" = paste( files.name, "- Results/CSV Sheets" ),
										"Summary Sheets Samples"= paste( files.name, "- Results/Summary Sheets - Samples" ),
										"Summary Sheets Integral"= paste( files.name, "- Results/Summary Sheets - Integral" )
	)
	
	
	dir.create( paste( files.name, "- Results" ), showWarnings=F )
	dir.create( paste( files.name, "- Results/CSV Sheets" ), showWarnings=F )

	

	if( !exists("colramp") )
		colramp = colorRampPalette( c( "darkblue","blue","lightblue","green","yellow","red","darkred" ) )
	
	
	
	
	if( exists("group.labels") && length(group.labels) != ncol(indata) || 
			exists("group.colors") && length(group.colors) != ncol(indata) )
	{
		rm(group.labels)
		rm(group.colors)
		cat( "\n!!Group assignment doesnt fit number of samples!!\n" ); flush.console()			
	}
	
	if( exists("group.labels") && max( table(group.labels) ) == 1 )
	{
		rm(group.labels)
		if( exists("group.colors") ) rm(group.colors)
		cat( "\n!!Each sample has an onw group!!\n" ); flush.console()			
	}
	
	if( exists("group.labels") )
	{	
		group.labels = as.character( group.labels )
		names(group.labels) = colnames( indata )
                                                       
		if( !exists("group.colors") )                 
		{
			group.colors = rep( "", ncol( indata ) )
			for( i in 1:length( unique( group.labels ) ) )
			{
				group.colors[ which(group.labels == unique( group.labels )[i] ) ] = colorRampPalette( c("blue3","blue","lightblue","green","gold","red","red3") )(length( unique( group.labels ) ))[i]				
			}
		}	
		
		if( length( unique( substr( group.colors, 1, 1 ) ) > 1 ) || unique( substr( group.colors, 1, 1 ) )[1] != "#" )	# catch userdefined group.colors --> convert to #hex
		{
			group.colors = apply( col2rgb( group.colors ), 2, function(x){ rgb(x[1]/255,x[2]/255,x[3]/255) } )
		}

		names(group.colors) = colnames( indata )
		
	} else
	{
		group.labels = rep( "sample", ncol( indata ) )
		names(group.labels) = colnames( indata )

		group.colors = colramp( ncol( indata )  )
		names(group.colors) = colnames( indata )
	}

	unique.group.colors = group.colors[ match( unique(group.labels), group.labels ) ]
	names(unique.group.colors) = unique(group.labels)

	


	LETTERS = c( LETTERS,   as.vector( sapply( 1:10, function(x){ paste( LETTERS, x, sep=""  ) } ) )   )
	letters = c( letters,   as.vector( sapply( 1:10, function(x){ paste( letters, x, sep=""  ) } ) )   )
	
	

	
	
	# prepare data

	
	indata.sample.mean = colMeans(indata)	
	
	source("lib/source/Quality Check.r")
	
	

	if( preferences$sample.quantile.normalization )
	{
		indata = Quantile.Normalization( indata )
		
		if( exists("indata.original") )
		{
			indata.original = Quantile.Normalization( indata.original )
			cat( "\n!!Separate quantile normalization of indata AND indata.original!!\n" ); flush.console()			
		}
	}	
	
	
	if( exists("indata.original") && any( dim(indata) != dim(indata.original) ) )
	{
		rm(indata.original)
		cat( "\n!!Existing 'indata.original' does not fit 'indata' object!!\n" ); flush.console()		
	}
	if( !exists("indata.original") )
	{
		indata.original = indata
	}
	
	
	
	
	if( preferences$error.model == "replicates" )
	{
		indata = do.call( cbind, by( t(indata), colnames(indata), mean ) )[ , unique( colnames(indata) ) ]

		group.labels = group.labels[ colnames(indata) ]
		group.colors = group.colors[ colnames(indata) ]
		
		indata.sample.mean = tapply( indata.sample.mean, colnames(indata.original), mean )[ colnames(indata) ]

	} else
	{
		colnames( indata ) = make.unique( colnames( indata ) )
		names( group.labels ) = make.unique( names( group.labels ) )
		names( group.colors ) = make.unique( names( group.colors ) )
	}





	indata.mean.level = rowMeans( indata )

	if( preferences$feature.mean.normalization )
	{		
		indata = indata - indata.mean.level
	}



	
	



	cat( "\n\nLoad Annotation Data\n\n" ); flush.console()
	source("lib/source/Prepare Annotation.r")


	
	

## SOM Processing


	cat( "\nProcessing SOM  " );	flush.console()


	
	som.result = som.init( indata, xdim=preferences$dim.som1, ydim=preferences$dim.som1, init="linear" )

	# Rotate/Flip First lvl SOMs

	if( preferences$rotate.som1 > 0 )
		for( i in 1:preferences$rotate.som1 )
		{
			o = matrix( c( 1:(preferences$dim.som1^2) ), preferences$dim.som1, preferences$dim.som1, byrow=T )
			o = o[ rev(1:preferences$dim.som1) ,]
			som.result = som.result[ as.vector(o), ]
		}
	if( preferences$flip.som1 )
		{
			o = matrix( c( 1:(preferences$dim.som1^2) ), preferences$dim.som1, preferences$dim.som1, byrow=T )
			som.result = som.result[ as.vector(o), ]
		}
		
	# Train SOM

	#format(Sys.time(), "%a %b %d %X\n" )

	t1 = system.time({
	som.result = som.train( indata, som.result, xdim=preferences$dim.som1, ydim=preferences$dim.som1, 
				alpha=0.05, radius=preferences$dim.som1, rlen=nrow(indata)*2*preferences$training.extension, inv.alp.c=nrow(indata)*2*preferences$training.extension/100 )
	})
	cat( "(remaining ~", ceiling(5*t1[3]/60),"min /",round(5*t1[3]/3600,1),"h )","\n" );	flush.console()
	
	som.result = som.train( indata, som.result$code, xdim=preferences$dim.som1, ydim=preferences$dim.som1, 
				alpha=0.02, radius=min(3, preferences$dim.som1), rlen=nrow(indata)*10*preferences$training.extension, inv.alp.c=nrow(indata)*10*preferences$training.extension/100 )

	#som.result = som( indata, xdim=preferences$dim.som1, ydim=preferences$dim.som1 )

	#format(Sys.time(), "%a %b %d %X\n" )


	metadata = som.result$code
	colnames(metadata) = colnames( indata )

	som.result$code = NA







	loglog.metadata = apply( metadata, 2, function(x)
	{
		meta.sign = sign( x )
		meta = log10( abs( x ) )
		meta = meta - min( meta, na.rm=T )
		return(meta * meta.sign)
	})


	WAD.metadata = apply( metadata ,2, function(x)
	{
		w = ( x - min( x ) ) / ( max( x ) - min( x ) )
		return(x * w)
	})


##	Group SOMs

	if( length(unique( group.labels ) ) > 1 )			# mean group metagenes
	{
		group.metadata = do.call( cbind, by( t(metadata), group.labels, colMeans ) )[,unique(group.labels)]

		loglog.group.metadata = do.call( cbind, by( t(loglog.metadata), group.labels, colMeans ) )[,unique(group.labels)]

		WAD.group.metadata = do.call( cbind, by( t(WAD.metadata), group.labels, colMeans ) )[,unique(group.labels)]
	}






## set up SOM dependent variables

	genes.coordinates = apply( som.result$visual[,c(1,2)]+1, 1, paste, collapse=" x " )
	names( genes.coordinates ) = rownames(indata)


	som.nodes = ( som.result$visual[,"x"] + 1 ) +  som.result$visual[,"y"] * preferences$dim.som1
	names( som.nodes ) = rownames( indata )












## Output




	save.image( paste( files.name, " pre.RData" , sep="" ) )



	cat( "Processing Differential Expression\n" ); flush.console()
	source("lib/source/Calc Statistics.r")


	cat( "Spot Detection\n" ); flush.console()
	source("lib/source/Detect Spots Samples.r")
	source("lib/source/Detect Spots Integral.r")


	source("lib/source/Group Assignment.r")


	cat( "Plotting Sample Portraits\n" ); flush.console()
	source("lib/source/Sample Expression Portraits.r")
	source("lib/source/Sample Rank Maps.r")


	cat( "Processing Supporting Information\n" ); flush.console()
	source("lib/source/Supporting Maps.r")
	source("lib/source/Entropy Profiles.r")	
	source("lib/source/Topology Profiles.r")


	cat( "Processing 2nd level Metagene Analysis\n" ); flush.console()
	dir.create( paste( files.name, "- Results/2nd lvl Metagene Analysis" ), showWarnings=F )
	
	source("lib/source/2nd lvl Similarity Analysis.r")
	source("lib/source/2nd lvl Correlation Analysis.r")
	source("lib/source/2nd lvl Component Analysis.r")
	source("lib/source/2nd lvl SOM.r")
	
	

	if( preferences$geneset.analysis )
	{
		dir.create( paste( files.name, "- Results/Geneset Analysis" ), showWarnings=F )
		
		source("lib/source/Geneset Statistic Samples.r")
		source("lib/source/Geneset Statistic Integral.r")
		source("lib/source/Geneset Overviews.r")
		source("lib/source/Geneset Profiles + Maps.r")

		source("lib/source/Cancer Hallmarks.r")
		source("lib/source/Chromosome Expression Reports.r")
	}



	cat( "Gene Lists\n" ); flush.console()
	source("lib/source/Gene Lists.r")

	cat( "Summary Sheets: Samples\n" ); flush.console()
	source("lib/source/Summary Sheets Samples.r")
	
	cat( "Summary Sheets: Spots\n" ); flush.console()
	source("lib/source/Summary Sheets Integral.r")
	
	

	cat( "Processing 3rd level Spot Analysis\n" ); flush.console()	
	dir.create( paste( files.name, "- Results/3rd lvl Spot Analysis" ), showWarnings=F )
	
	source("lib/source/3rd lvl Chromosomal Enrichment.r")	
	source("lib/source/3rd lvl Summary Sheets.r")		
	source("lib/source/3rd lvl Networks.r")	




	cat( "Generating HTML Report\n" ); flush.console()
	source("lib/source/HTML Summary.r")
	source("lib/source/HTML Sample Summary.r")
	source("lib/source/HTML Integral Summary.r")	
	source("lib/source/HTML Geneset Analysis.r")


	
	cat( "Clean and store Workspace\n" ); flush.console()	
	source("lib/source/Workspace Cleanup.r")
	save.image( paste( files.name, ".RData" , sep="" ) )
	
	
	if( file.exists( paste( files.name, " pre.RData" , sep="") ) && file.exists( paste( files.name, ".RData" , sep="") ) )		
		r = file.remove( paste( files.name, " pre.RData" , sep="") )


	
	
	## additional scripts

	source("lib/source/Group Analyses.r")

	source("lib/source/Difference Analyses.r")
	
	
	
#	cat( "Spot Filtering\n" ); flush.console()
#	source("lib/source/3rd lvl Overexpression Genenet.r")	
#	source("lib/source/3rd lvl Spot Filter.r")
	source("lib/source/Signature Sets.r")
	
	
	
	cat( "Finished:", format(Sys.time(), "%a %b %d %X\n\n" ) )
	flush.console()





	
	

