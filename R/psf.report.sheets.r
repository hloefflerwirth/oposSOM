
circle <- function(x, y, r, ...) {
  xx <- c()
  yy <- c()
  for(i in 1:100 ) {
    xx[i] <- cos(i*6.28/100)
    yy[i] <- sin(i*6.28/100)
  }
  polygon( r * xx + x, r * yy + y, ...)
}


plot.psf.titlepage <- function(env, psf.object, signal.values, bar.colors)
{
  layout(matrix(c(1,2,3,4,4,4,5,5,5),3,byrow=TRUE))
  
  par(mar=c(0,0,0,0))
  plot(0, type="n", axes=FALSE, xlab="", ylab="", xlim=c(0,1), ylim=c(0,1))
  text(0.05, 0.94, psf.object$attrs$title , cex=2, adj=0)
  
  ### Population map of all nodes ###
  n.map <- matrix(0,env$preferences$dim.1stLvlSom,env$preferences$dim.1stLvlSom)
  pw.genes <- unlist(sapply(psf.object$graph@nodeData@data,function(x)x$genes))
  pw.genes <- unique(env$gene.info$ensembl.mapping[which(env$gene.info$ensembl.mapping$ensembl_gene_id %in% pw.genes),1] )
  pw.metagenes <-  env$som.result$feature.BMU[pw.genes]
  n.map[as.numeric(names(table(pw.metagenes)))] <- table(pw.metagenes)
  n.map[which(n.map==0)] <- NA
  n.map <- matrix(n.map, env$preferences$dim.1stLvlSom)
  
  lim <- c(1,env$preferences$dim.1stLvlSom) + env$preferences$dim.1stLvlSom * 0.01 * c(-1, 1)
  colr <- env$color.palette.portraits(1000)[(na.omit(as.vector(n.map)) - min(n.map,na.rm=TRUE)) /
                                          max(1, (max(n.map,na.rm=TRUE) - min(n.map,na.rm=TRUE))) *
                                          999 + 1]
  
  par(mar=c(3,6,3,6))
  spot.background <- get(paste("spot.list.",env$preferences$standard.spot.modules,sep=""),envir = env)$overview.mask
  image(matrix(spot.background, env$preferences$dim.1stLvlSom), col="gray90", axes=FALSE )
  par(new=TRUE)
  plot(which(!is.na(n.map), arr.ind=TRUE), xlim=lim, ylim=lim, pch=16, axes=FALSE,
       xlab="",ylab="", xaxs="i", yaxs="i", col=colr, main="all genes",
       cex=0.5 + na.omit(as.vector(n.map)) / max(n.map,na.rm=TRUE) * 2.8)
  title(sub=paste("maximum =", max(n.map,na.rm=TRUE)),line=0)
  box()
  
  ### Population map of sink-nodes ###
  if(length(psf.object$sink.nodes)>0)
  {
    n.map <- matrix(0,env$preferences$dim.1stLvlSom,env$preferences$dim.1stLvlSom)
    pw.genes <- unlist(sapply(psf.object$graph@nodeData@data[psf.object$sink.nodes],function(x)x$genes))
    pw.genes <- unique(env$gene.info$ensembl.mapping[which(env$gene.info$ensembl.mapping$ensembl_gene_id %in% pw.genes),1] )
    pw.metagenes <- env$som.result$feature.BMU[pw.genes]
    if(length(pw.metagenes)>0)
    {
      n.map[as.numeric(names(table(pw.metagenes)))] <- table(pw.metagenes)
      n.map[which(n.map==0)] <- NA
      n.map <- matrix(n.map, env$preferences$dim.1stLvlSom)
      
      lim <- c(1,env$preferences$dim.1stLvlSom) + env$preferences$dim.1stLvlSom * 0.01 * c(-1, 1)
      colr <- env$color.palette.portraits(1000)[(na.omit(as.vector(n.map)) - min(n.map,na.rm=TRUE)) /
                                              max(1, (max(n.map,na.rm=TRUE) - min(n.map,na.rm=TRUE))) *
                                              999 + 1]
      
      par(mar=c(3,6,3,6))
      
      image(matrix(spot.background, env$preferences$dim.1stLvlSom), col="gray90", axes=FALSE )
      par(new=TRUE)
      plot(which(!is.na(n.map), arr.ind=TRUE), xlim=lim, ylim=lim, pch=16, axes=FALSE,
           xlab="",ylab="", xaxs="i", yaxs="i", col=colr, main="sink node genes",
           cex=0.5 + na.omit(as.vector(n.map)) / max(n.map,na.rm=TRUE) * 2.8)
      title(sub=paste("maximum =", max(n.map,na.rm=TRUE)),line=0)
      box()
      
    } else frame()
    
    
    ### Profile of mean sink signals ###
    ylim <- c(-2, 2)
    par(mar=c(2,7,4,5))
    
    barplot( sapply(signal.values,function(x) mean(log10(x$signal.at.sinks),na.rm=T) ),
             beside=TRUE, col=bar.colors, names.arg=rep("",length(signal.values)),
             ylim=ylim, border=if (length(signal.values) < 80) "black" else NA )
    mtext(bquote("<log"[10] ~ "s>"), side=2, line=2.5, cex=1.5)
    
    ### Profile of max sink signals ###
    ylim <- c(-2, 2)
    par(mar=c(6,7,0,5))
    
    bar.coords <- barplot( sapply(signal.values,function(x) max(log10(x$signal.at.sinks),na.rm=T) ),
                           beside=TRUE, names.arg=rep("",length(signal.values)),
                           col=bar.colors, ylim=ylim, border=if (length(signal.values) < 80) "black" else NA )
    mtext(bquote("log"[10] ~ "s"[max]), side=2, line=2.5, cex=1.5)
    
    if (length(signal.values)<100)
      text(bar.coords, par('usr')[3], labels=names(signal.values), srt=45, adj=c(1.1,1.1), xpd=TRUE)
  }
}


plot.psf.pathway.keggrest <- function( kegg.pathway, signal.values=NULL, signal.values.lim, main="", highlight.genes=NULL, color.palette=NULL )
{
	width = ncol(kegg.pathway$pathway.img)
	height = nrow(kegg.pathway$pathway.img)
	max.xy = max(width,height)

	par(mar = c(0, 0, 0, 0))
	plot(c(0-(max.xy-width)/2, width+(max.xy-width)/2), c(0-(max.xy-height)/2, height+(max.xy-height)/2), type = "n", xlab = "", ylab = "", xaxs = "i", yaxs = "i", axes=F)

	rasterImage(kegg.pathway$pathway.img, 0, 0, width, height, interpolate = F)

	title.node <- which( sapply( kegg.pathway$pathway.info, function(x) grepl("TITLE:", x$graphics$label.long ) ) )
	if(length(title.node)>0)
	{
		rect( kegg.pathway$pathway.info[[title.node]]$graphics$x-kegg.pathway$pathway.info[[title.node]]$graphics$width*0.6,
					height-kegg.pathway$pathway.info[[title.node]]$graphics$y+kegg.pathway$pathway.info[[title.node]]$graphics$height*0.6,
					kegg.pathway$pathway.info[[title.node]]$graphics$x+kegg.pathway$pathway.info[[title.node]]$graphics$width*0.6,
					height-kegg.pathway$pathway.info[[title.node]]$graphics$y-kegg.pathway$pathway.info[[title.node]]$graphics$height*0.6,
					col="gray90", lwd=0.5 )

		text( kegg.pathway$pathway.info[[title.node]]$graphics$x, height - kegg.pathway$pathway.info[[title.node]]$graphics$y, main,
					cex = .5, col = "black", family="sans" )
	}
	
	
	for( i in names( which( sapply(kegg.pathway$pathway.info,"[[", "type" ) == "gene" ) ) )
	{
		node.col = "gray90"
		if( paste( kegg.pathway$pathway.info[[i]]$id, kegg.pathway$pathway.info[[i]]$graphics$x, kegg.pathway$pathway.info[[i]]$graphics$y ) %in% highlight.genes )
		{
			node.col = "olivedrab1"

		} else if( i %in% names(signal.values) )
		{
			node.col = color.palette(1000)[999*(signal.values[i]-signal.values.lim[1])/(signal.values.lim[2]-signal.values.lim[1])+1]
		}

		rect( kegg.pathway$pathway.info[[i]]$graphics$x-kegg.pathway$pathway.info[[i]]$graphics$width*0.5,
					height-kegg.pathway$pathway.info[[i]]$graphics$y+kegg.pathway$pathway.info[[i]]$graphics$height*0.5,
					kegg.pathway$pathway.info[[i]]$graphics$x+kegg.pathway$pathway.info[[i]]$graphics$width*0.5,
					height-kegg.pathway$pathway.info[[i]]$graphics$y-kegg.pathway$pathway.info[[i]]$graphics$height*0.5,
					col=node.col, lwd=0.5 )

		text( kegg.pathway$pathway.info[[i]]$graphics$x, height - kegg.pathway$pathway.info[[i]]$graphics$y, kegg.pathway$pathway.info[[i]]$graphics$label.short,
					cex = ifelse( nchar(kegg.pathway$pathway.info[[i]]$graphics$label.short)<=6,.3,.24), col = "black", family="sans" )
	}
	
	for( i in which( sapply(kegg.pathway$pathway.info,"[[", "type" ) == "gene" ) )
	  if( i %in% names(signal.values) )	  
  	{
  	  node.col = color.palette(1000)[999*(signal.values[i]-signal.values.lim[1])/(signal.values.lim[2]-signal.values.lim[1])+1]
  
  	  rect( kegg.pathway$pathway.info[[i]]$graphics$x-kegg.pathway$pathway.info[[i]]$graphics$width*0.5,
  	        height-kegg.pathway$pathway.info[[i]]$graphics$y+kegg.pathway$pathway.info[[i]]$graphics$height*0.5,
  	        kegg.pathway$pathway.info[[i]]$graphics$x+kegg.pathway$pathway.info[[i]]$graphics$width*0.5+2,
  	        height-kegg.pathway$pathway.info[[i]]$graphics$y-kegg.pathway$pathway.info[[i]]$graphics$height*0.5,
  	        col=node.col, lwd=0.5 )
  
  	  text( kegg.pathway$pathway.info[[i]]$graphics$x, height - kegg.pathway$pathway.info[[i]]$graphics$y, kegg.pathway$pathway.info[[i]]$graphics$label.short,
  	        cex = 1-0.02*nchar(kegg.pathway$pathway.info[[i]]$graphics$label.short), col = "black", family="sans" )
  	}

	for( i in which( sapply(kegg.pathway$pathway.info,"[[", "type" ) == "compound" ) )
	  if( i %in% names(signal.values) )	  
  	{
  	  node.col = color.palette(1000)[999*(signal.values[i]-signal.values.lim[1])/(signal.values.lim[2]-signal.values.lim[1])+1]
  
  	  circle( kegg.pathway$pathway.info[[i]]$graphics$x, height-kegg.pathway$pathway.info[[i]]$graphics$y,
  	          kegg.pathway$pathway.info[[i]]$graphics$width*0.75,
  	          col=node.col, lwd=0.5 )
  	}

	map.nodes <- which( sapply(kegg.pathway$pathway.info,"[[", "type" ) == "map" )
	map.nodes <- setdiff( map.nodes, title.node )

	for( i in map.nodes )
	{
	  rect( kegg.pathway$pathway.info[[i]]$graphics$x-kegg.pathway$pathway.info[[i]]$graphics$width*0.5-1,
	        height-kegg.pathway$pathway.info[[i]]$graphics$y+kegg.pathway$pathway.info[[i]]$graphics$height*0.5+1,
	        kegg.pathway$pathway.info[[i]]$graphics$x+kegg.pathway$pathway.info[[i]]$graphics$width*0.5+2,
	        height-kegg.pathway$pathway.info[[i]]$graphics$y-kegg.pathway$pathway.info[[i]]$graphics$height*0.5-3,
	        col = "gray95", lwd=0.5 )
	  
	  text( kegg.pathway$pathway.info[[i]]$graphics$x, height - kegg.pathway$pathway.info[[i]]$graphics$y, kegg.pathway$pathway.info[[i]]$graphics$label.short,
	        cex = 1-0.02*nchar(kegg.pathway$pathway.info[[i]]$graphics$label.short), col = "black", family="sans" )
	}

}



psf.report.sheets <- function(env, psf.results, output.path, bar.colors )
{
  util.info("Writing: PSF report sheets")
  progressbar <- newProgressBar(min = 0, max = length(kegg.collection)); cat("\r")
  
  for(pw in 1:length(kegg.collection) )
  {

    kegg.collection[[pw]]$id = strsplit( kegg.collection[[pw]]$attrs$name, ":" )[[1]][2]
    data(list=kegg.collection[[pw]]$id)
    

    filename <- file.path(output.path, paste( make.names(names(kegg.collection)[pw]),".pdf",sep="") )
    pdf(filename, 29.7/2.54, 21/2.54)

    plot.psf.titlepage( env, psf.object=kegg.collection[[pw]], signal.values=psf.results[[pw]], bar.colors )
                        

    par(mfrow=c(1,1),mar=c(0,0,0,0))

    node.genes <- lapply(kegg.collection[[pw]]$graph@nodeData@data,function(x)x$genes)
    node.genes <- lapply(node.genes,function(x) unique(env$gene.info$ensembl.mapping[which(env$gene.info$ensembl.mapping$ensembl_gene_id %in% x),1] ) )
    node.genes <- sapply(node.genes,function(x) any(!is.na(x)) )
    node.genes <- sapply( kegg.collection[[pw]]$graph@nodeData@data[which(node.genes)], function(x) paste(x$kegg.name, x$kegg.gr.x, x$kegg.gr.y ) )
 
    plot.psf.pathway.keggrest( kegg.pathway = kegg.data, highlight.genes=node.genes, main=paste(names(kegg.collection)[pw],"\ngenes with data" ) )
 

    par(mfrow=c(1,1),mar=c(0,0,0,0))
    
    sink.genes <- sapply(kegg.collection[[pw]]$graph@nodeData@data,function(x)x$kegg.id %in% kegg.collection[[pw]]$sink.nodes )
    sink.genes <- sapply( kegg.collection[[pw]]$graph@nodeData@data[which(sink.genes)], function(x) paste(x$kegg.name, x$kegg.gr.x, x$kegg.gr.y ) )
    
    plot.psf.pathway.keggrest( kegg.pathway = kegg.data, highlight.genes=sink.genes, main=paste(names(kegg.collection)[pw],"\nsink nodes") )
 
 
    for( m in 1:length(psf.results[[pw]]) )
    {
      signal.values = log10( psf.results[[pw]][[m]]$signal.at.nodes )
      names(signal.values) = sapply( kegg.collection[[pw]]$graph@nodeData@data, function(x) paste(x$kegg.name, x$kegg.gr.x, x$kegg.gr.y ) )

      plot.psf.pathway.keggrest( kegg.pathway = kegg.data, signal.values = signal.values,
                                  signal.values.lim = c(-1,1)*max( abs( log10( sapply( psf.results[[pw]], function(x) x$signal.at.nodes ) ) ) ),
                                  main=paste(names(kegg.collection)[pw],"\n",names(psf.results[[pw]])[m] ),
                                  color.palette=env$color.palette.heatmaps)
    }

    dev.off()

    setTxtProgressBar( progressbar, progressbar$getVal()+1 )
  }
  progressbar$kill()    

}

