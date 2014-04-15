
	dir.create( paste( files.name, "- Results/CSV Sheets/Chromosomal Enrichment" ), showWarnings=F )


	sort.label = function( x )
	{
		numeric = function(x) suppressWarnings( as.numeric(x) )		
		nonnumeric = function(x)	ifelse(is.na(suppressWarnings(as.numeric(x))), toupper(x), NA)		
		
		x <- as.character(x)	
		delim="\\$\\@\\$"	
		delimited <- gsub("([+-]{0,1}[0-9]+([eE][\\+\\-]{0,1}[0-9]+){0,1})",	paste(delim,"\\1",delim,sep=""), x)
		
		step1 <- strsplit(delimited, delim)
		step1 <- lapply( step1, function(x) x[x>""] )
		
		step1.numeric <- lapply( step1, numeric )
		step1.character <- lapply( step1, nonnumeric )
		
		maxelem <- max(sapply(step1, length))
		
		step1.numeric.t <- lapply( 1:maxelem, function(i)	sapply(step1.numeric, function(x)x[i]) )
		step1.character.t <- lapply( 1:maxelem, function(i)	sapply(step1.character, function(x)x[i]) ) 
		
		rank.numeric   <- sapply(step1.numeric.t,rank)
		rank.character <- sapply(step1.character.t, function(x) as.numeric(factor(x)))
		rank.numeric[!is.na(rank.character)] <- 0	
		rank.character <- t( t(rank.character) + apply(matrix(rank.numeric),2,max,na.rm=TRUE) )
		rank.overall <- ifelse(is.na(rank.character),rank.numeric,rank.character)
		
		order.frame <- as.data.frame(rank.overall)
		x = x[ do.call("order",order.frame) ]
		return( x )
	}

	



	
	plot.set.list.chromosomes = function( set.list, main )
	{
	
		sorted.unique.gene.positions = sort.label( unique(gene.positions) )

		n.all.chr.genes = length(unlist(gene.positions.list))
		
		
		intersect.counts = matrix( 0, length(unique(gene.positions)), length(set.list$spots), dimnames=list( sorted.unique.gene.positions, names(set.list$spots) ) )
	
		for( m in names( set.list$spots ) )
		{
			intersect.counts.spot = sapply( gene.positions.list, function(x)
			{
				sapply( x, function(y)
				{				
#					-log10( fisher.test( matrix(c(n.all.chr.genes-length(set.list$spots[[ m ]]$genes)-length(y)-length( intersect( set.list$spots[[ m ]]$genes,  y ) ),length(set.list$spots[[ m ]]$genes),length(y),length( intersect( set.list$spots[[ m ]]$genes,  y ) )),2) )$p.value )
					length( intersect( set.list$spots[[ m ]]$genes,  y ) ) / length(y)
				} )
			} )
			
			intersect.counts.spot = unlist( intersect.counts.spot )
			names(intersect.counts.spot) = sub( ".", " ", names(intersect.counts.spot), fixed=T )
						
			intersect.counts[ names(intersect.counts.spot), m ] = intersect.counts.spot
		}
		
		write.csv2( intersect.counts, file=paste( files.name, " - Results/CSV Sheets/Chromosomal Enrichment/", main, ".csv" ,sep="") )
		

		
		colkey= c( colorRampPalette( c( rep("white",2),"green",rep("yellow2",1),rep("orange",2),rep("red",2),rep("darkred",2),rep("purple",1) ) )(2000),
							 colorRampPalette( c("purple","black") )(3000), rep("black",5000 ) )
		

		
		chr = sapply( strsplit( rownames(intersect.counts), " ", fixed=T ), head, 1 )
		chr.sep.level = table( chr )[unique(chr)]
		chr.sep.level = rev(chr.sep.level)	
		chr.sep.level = chr.sep.level[-length(chr.sep.level)]
		chr.sep.level = sapply( 1:length(chr.sep.level), function(i) sum(chr.sep.level[1:i]) )		

		chr.lab.level = table( chr )[unique(chr)]
		chr.lab.level = rev(chr.lab.level)	
		chr.lab.level = c( 0, sapply( 1:(length(chr.lab.level)-1), function(i) sum(chr.lab.level[1:i]) ) ) + chr.lab.level/2
		
		
		par(mfrow=c(1,1), mar=c(2,3,2,2) )
		image( 1:ncol(intersect.counts), 1:nrow(intersect.counts), t(intersect.counts[nrow(intersect.counts):1,]), axes=F, xlab="", ylab="", col=colkey, main=main, zlim=c(0,1) )
			box()
			abline( h=chr.sep.level+0.5, lwd=0.1, lty=3 )
		
			axis( 1, 1:ncol(intersect.counts), labels = colnames(intersect.counts), las = 1, line = -0.5, tick = 0, cex.axis=0.8 )
			axis( 2, 1:nrow(intersect.counts), labels = rev(rownames(intersect.counts)), las = 2, line = -0.5, tick = 0, cex.axis=0.24 )
			axis( 2, chr.lab.level,	rev(unique(chr)), las=2, tick=F, line=0.6 )
		
		
		par(new=T, mar=c(1,0,0,0) );	layout( matrix( c(0,0,0,0,1,0,0,0,0), 3, 3 ), c( 1.6, 0.03, 0.005 ), c( 0.115, 1.15, 2.15 ) )
		image( matrix( 1:1000, 1, 1000 ), col = colkey, axes=F )
			axis( 2, at=seq( 0, 1, 0.1 ), paste( seq( 0, 100, 10 ),"%",sep=""), las=2, tick=F, pos=0.5, cex.axis=0.6 )		
			box()
		


		
		source("lib/source/3rd lvl Chromosome Association Maps.r",local=T)
		
		
		
		for( m in 1:length(set.list$spots) )
		{

			intersect.counts = sapply( gene.positions.list, function(x)
			{
				sapply( x, function(y)
				{
					length( intersect( set.list$spots[[ m ]]$genes,  y ) ) / length(y)
				} )
			} )
			names(intersect.counts) = names(gene.positions.list)
		
			l = c( rep(1,ceiling( length(gene.positions.list)/2 )), 2:(length(gene.positions.list)+1) )
			l = c( l, rep( 0, ( ceiling( length(gene.positions.list)/2 )*3 ) - length(l) ) )
			layout( matrix( l, 3, ceiling( length(gene.positions.list)/2 ), byrow=T ), heights=c(0.05,1,1) )
			
			par( mar=c(0,0,0,0) )
			plot( 0, type="n", xlab="", ylab="", axes=F, xlim=c(0,1) )
			text( 0.005, 0, paste("Spot",names(set.list$spots)[m]), cex=2 )			
			
			
			par( mar=c(1.5,2,0.2,0.3) )
			
			for( chromosome in sort.label(names(gene.positions.list)) )
			{
				x = intersect.counts[[chromosome]]
				
#				image( 1, 1:length(x),  matrix( rev(x), 1, length(x) ), zlim=c( 0, max( unlist( intersect.counts ) ) ), axes=F, xlab="", ylab="", col=colorRampPalette( c("white","yellow","red") )(1000) )
				image( 1, 1:length(x),  matrix( rev(x), 1, length(x) ), zlim=c( 0, 1), axes=F, xlab="", ylab="", col=colkey )				
				box()
				axis(1, 1, labels = chromosome, las = 1, line = -0.8, tick = 0, cex.axis=1.2 )
				axis(2, 1:length(x), labels = rev(names(x)), las = 2, line = -0.5, tick = 0, cex.axis=1 )
			}
						
		}

	}




	pdf( paste( files.name, " - Results/3rd lvl Spot Analysis/Overexpression Chromosomal Enrichment.pdf", sep="" ), 21/2.54, 29.7/2.54 )

	if( length( gene.positions.list ) > 0 )
		plot.set.list.chromosomes( set.list=GS.infos.overexpression, main="Overexpression Spot Chromosome Map" )

	dev.off()


	pdf( paste( files.name, " - Results/3rd lvl Spot Analysis/Underexpression Chromosomal Enrichment.pdf", sep="" ), 21/2.54, 29.7/2.54 )

	if( length( gene.positions.list ) > 0 )
		plot.set.list.chromosomes( GS.infos.underexpression, main="Underexpression Spot Chromosome Map" )

	dev.off()


	pdf( paste( files.name, " - Results/3rd lvl Spot Analysis/K-Means Cluster Chromosomal Enrichment.pdf", sep="" ), 21/2.54, 29.7/2.54 )
	
	if( length( gene.positions.list ) > 0 )
		plot.set.list.chromosomes( GS.infos.kmeans, main="K-Means Cluster Chromosome Map" )
	
	dev.off()

	pdf( paste( files.name, " - Results/3rd lvl Spot Analysis/Group Overexpression Chromosomal Enrichment.pdf", sep="" ), 21/2.54, 29.7/2.54 )
	
	if( length( gene.positions.list ) > 0 )
		plot.set.list.chromosomes( set.list=GS.infos.group.overexpression, main="Group Overexpression Chromosome Map" )
	
	dev.off()








