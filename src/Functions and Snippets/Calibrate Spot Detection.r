
###########################################################################

#################### 			overexpression spots 		 ########################

###########################################################################

# 	# Fine tune parameter old "flooding hills spot detection"
# 	preferences$summary.spot.minsize
# 	preferences$summary.spot.maxsize
# 	preferences$summary.spot.startlevel
# 
# 	# preferences$summary.spot.minsize = 4
# 	# preferences$summary.spot.maxsize = 10
# 	# preferences$summary.spot.startlevel = 0.9
# 
# 	preferences$summary.spot.minsize = 1
# 	preferences$summary.spot.maxsize = 4
# 	preferences$summary.spot.startlevel = 0.9
# 
# 
# 	source("Functions and Snippets/f - Find  Overexpression spots - old.r")



	preferences$summary.spot.threshold
	preferences$summary.spot.core
	preferences$group.spot.core
	preferences$group.spot.threshold

	# preferences$summary.spot.threshold = 0.95
	# preferences$summary.spot.core = 5
	# preferences$group.spot.core = 5
	# preferences$group.spot.threshold = 0.75

	preferences$summary.spot.threshold = 0.98
	preferences$summary.spot.core = 6
	preferences$group.spot.core = 3
	preferences$group.spot.threshold = 0.88

	source("lib/f - Detect Spots Integral.r")
	
	colramp = colorRampPalette( c( "darkblue","blue","lightblue","green","yellow","red","darkred" ) )

	
	par(mfrow=c(2,2), mar=c(2,2,2,2))
	image( matrix( GS.infos.overexpression$overview.map, preferences$dim.som1 ), col=colramp(1000) )
	image( matrix( GS.infos.overexpression$overview.mask, preferences$dim.som1 ), col=colramp(1000) )
	par( new=T )
	plot( 0, type="n", axes=F, xlab="", ylab="", xlim=c(0,preferences$dim.som1), ylim=c(0,preferences$dim.som1), xaxs="i", yaxs="i" )
	points( do.call( rbind, lapply( GS.infos.overexpression$spots, function(x) x$position ) ), pch=16, cex=1, col="black" )		
	points( do.call( rbind, lapply( GS.infos.overexpression$spots, function(x) x$position ) ), pch=1, cex=1, col="white" )			
	
	length(GS.infos.overexpression$spots)


	image( matrix( GS.infos.group.overexpression$overview.map, preferences$dim.som1 ), col=colramp(1000) )
	image( matrix( GS.infos.group.overexpression$overview.mask, preferences$dim.som1 ), col=colramp(1000) )
	par( new=T )
	plot( 0, type="n", axes=F, xlab="", ylab="", xlim=c(0,preferences$dim.som1), ylim=c(0,preferences$dim.som1), xaxs="i", yaxs="i" )
	points( do.call( rbind, lapply( GS.infos.group.overexpression$spots, function(x) x$position ) ), pch=16, cex=1, col="black" )		
	points( do.call( rbind, lapply( GS.infos.group.overexpression$spots, function(x) x$position ) ), pch=1, cex=1, col="white" )			
	
	length(GS.infos.group.overexpression$spots)
	
	
	




# Re-run necessary parts of the pipeline

source("lib/f - Import.r")


dir.create( paste( files.name, "- Results" ), showWarnings=F )
dir.create( paste( files.name, "- Results/CSV Sheets" ), showWarnings=F )


if( preferences$geneset.analysis )
{
	dir.create( paste( files.name, "- Results/Geneset Analysis" ), showWarnings=F )
	source("lib/f - Geneset Statistic Integral.r")
	source("lib/f - Geneset Overviews.r")
	source("lib/f - Geneset Profiles + Maps.r")
}

source("lib/f - Gene Lists.r")
source("lib/f - Summary Sheets Integral.r")

dir.create( paste( files.name, "- Results/3rd lvl Spot Analysis" ), showWarnings=F )
source("lib/f - 3rd lvl Chromosomal Enrichment.r")	
source("lib/f - 3rd lvl Summary Sheets.r")		
source("lib/f - 3rd lvl Networks.r")	

source("lib/f - HTML Integral Summary.r")	

source("lib/f - Workspace Cleanup.r")
save.image( paste( files.name, ".RData" , sep="" ) )


#source("lib/f - 3rd lvl Overexpression Genenet.r")	
source("lib/f - Signature Sets.r")



	


