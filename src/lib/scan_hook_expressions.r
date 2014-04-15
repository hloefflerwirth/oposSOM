#suppressWarnings(


load.measures = function( file ) 
{
	suppressWarnings( { ret = read.table( file, sep="\t", skip=1 ) } )

	n = ret[ , 1 ]
	ret = as.numeric( as.character( ret[ 1:nrow(ret), 2 ] ) )

	names( ret ) = n[ 1:length(ret) ]

	return( ret )
}


load.variance = function( file ) 
{
	ret = read.table( file, sep="\t", skip=1 )

	n = ret[ , 1 ]
	ret = as.numeric( as.character( ret[ , 6 ] ) ) ^ 2

	names( ret ) = n

	return(  ret )
}



get.chiptype = function( folder, hook.version )
{
	ret = "NA"

	if( hook.version == "V1" )
	{
		ret = strsplit( readLines( paste( folder,"/statistics.dat",sep="") )[23], " = ", fixed=T )[[1]][2]
	}
	if( hook.version == "V1b" )
	{
		ret = strsplit( readLines( paste( folder,"/statistics.dat",sep="") )[20], " = ", fixed=T )[[1]][2]
	}
	if( hook.version == "V2" )
	{
		typesplit = strsplit( readLines( paste( folder,"/parameters.log",sep="") )[4], "/", fixed=T )[[1]]
		typesplit = typesplit[ length(typesplit) ]
		ret = sub( "_probe_tab", "", typesplit )
	}

	return( ret )
}



get.perc.N = function( folder, hook.version )
{
	if( hook.version == "V1" )
	{
		ret = NA
	}
	if( hook.version == "V1b" )
	{
		ret = list( perc.N = strsplit( readLines( paste( folder,"/statistics.dat",sep="") )[12], " = ", fixed=T )[[1]][2] )
	}
	if( hook.version == "V2" )
	{
		ret = NA
	}

	return( as.numeric( ret ) )
}




Scan.Hook.Expressions = function( path, selection )
{

	data.list = list()
	variance.list = list()
	hook.versions = c()
	chip.types = c()
	perc.N = c()

	l = list.files( as.character(path) )

	if( !missing( selection ) ) l = intersect( selection, l )

	if( length(l) == 0 ) stop("Wrong or empty path!")	


	cat( "Read Measures: ", path, "\n|" )
	for( i in 1:50 ) cat(" ")
	cat("|\n|")
	flush.console()

                                    
	for( i in 1:length(l) )
	{
		if ( file.exists( paste( path,"/",l[i],"/","ExpressionMeasures.txt",sep="") ) )
		{
			data.list[[ l[i] ]] = log10( load.measures( paste( path,"/",l[i],"/","ExpressionMeasures.txt",sep="") ) )

			hook.versions[ l[i] ] = "V1"
			chip.types[ l[i] ] = get.chiptype( paste( path,"/",l[i],sep=""), hook.versions[ l[i] ] )
			perc.N[ l[i] ] = get.perc.N( paste( path,"/",l[i],sep=""), hook.versions[ l[i] ] )
		}

		if ( file.exists( paste( path,"/",l[i],"/","Expression measures (gLog).res",sep="") ) )
		{
			data.list[[ l[i] ]] = load.measures( paste( path,"/",l[i],"/","Expression measures (gLog).res",sep="") )

			hook.versions[ l[i] ] = "V1b"
			chip.types[ l[i] ] = get.chiptype( paste( path,"/",l[i],sep=""), hook.versions[ l[i] ] )
			perc.N[ l[i] ] = get.perc.N( paste( path,"/",l[i],sep=""), hook.versions[ l[i] ] )
		}

		if ( file.exists( paste( path,"/",l[i],"/","ExpressionMeasure-PmExp10.dat",sep="") ) )
		{
			data.list[[ l[i] ]] = log10( load.measures( paste( path,"/",l[i],"/","ExpressionMeasure-PmExp10.dat",sep="") ) )
			if ( file.exists( paste( path,"/",l[i],"/","probesetVariance",sep="") ) )
				variance.list[[ l[i] ]] = load.variance( paste( path,"/",l[i],"/","probesetVariance",sep="") )

			hook.versions[ l[i] ] = "V2"
			chip.types[ l[i] ] = get.chiptype( paste( path,"/",l[i],sep=""), hook.versions[ l[i] ] )
			perc.N[ l[i] ] = get.perc.N( paste( path,"/",l[i],sep=""), hook.versions[ l[i] ] )
		}

		
		out.intervals = round( seq( 1, length(l), length(l) * 0.02 ) )[-1]
		out.intervals[ length( out.intervals ) : 49 ] = out.intervals[ length( out.intervals ) ]
		cat( paste( rep("#",length( which( out.intervals == i ) ) ), collapse="" ) );	flush.console()

	}


	min.length = min( sapply( data.list, length ) )

	data.list = lapply( data.list, function(x){ x = x[1:min.length] } )
	variance.list = lapply( variance.list, function(x){ x = x[1:min.length] } )



	data = matrix( unlist( data.list ), min.length, length(data.list) )
	colnames( data ) = names( data.list ) 
	rownames( data ) = names( data.list[[1]] )


	variance = NA
	if( length(variance.list)>0 )
	{
		variance = matrix( unlist( variance.list ), min.length, length(variance.list) )
		colnames( variance ) = names( variance.list ) 
		rownames( variance ) = names( variance.list[[1]] )
	}


	cat("#|\n")
	flush.console()


	return( list( Data.Matrix=data, Variance.Matrix=variance, Hook.Versions=hook.versions, Chip.Types=chip.types, Percent.N=perc.N  ) )
}
