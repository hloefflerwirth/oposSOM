

### Overexpression spots


thresh.global = max( max(mean.sp.gr.FC), -min(mean.sp.gr.FC) ) / 6


	sample.number = apply( mean.sp.FC, 1, function(x)
	{
		sum( x > thresh.global ) + sum( x < -thresh.global )
	} )
	group.number = apply( mean.sp.gr.FC, 1, function(x)
	{
		sum( x > thresh.global ) + sum( x < -thresh.global )
	} )
	

	
	
	if( any( group.number > 0 ) )
	{
	
		GS.infos.overexpression$filtered = T
	#	GS.infos.overexpression$spots = GS.infos.overexpression$spots[ which( sample.number >= 21 ) ]  # min(table(group.labels))*0.5
		GS.infos.overexpression$spots = GS.infos.overexpression$spots[ which( group.number > 0 ) ]
		GS.infos.overexpression$overview.mask[ -unique( unlist( sapply( GS.infos.overexpression$spots, function(x) x$metagenes ) ) ) ] = NA	
	
		
		source("lib/source/3rd lvl Overexpression Summary.r")
		source("lib/source/3rd lvl Overexpression Networks.r")	
		source("lib/source/3rd lvl Overexpression Genenet.r")

	
		pdf( paste( files.name, " - Results/Summary Sheets - Integral/Overexpression filtered.pdf", sep="" ), 29.7/2.54, 21/2.54 )
		plot.set.list( set.list=GS.infos.overexpression, main="Sample-Overexpression" )	
		dev.off()	
	 	
	}

