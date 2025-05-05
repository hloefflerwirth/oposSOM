
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
  text(0.05, 0.94, psf.object$title, cex=2, adj=0)
  
  ### Population map of all nodes ###
  n.map <- matrix(0,env$preferences$dim.1stLvlSom,env$preferences$dim.1stLvlSom)
  pw.genes <- unlist( sapply( psf.object$node.info, "[[", "ens.ids" ) ) 
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
    pw.genes <- unlist(sapply(psf.object$node.info[psf.object$sink.nodes],"[[","ens.ids"))
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
    
    ylim <- range( log10(signal.values) )
    
    par(mar=c(2,7,4,5))
    
    barplot( log10( colMeans(signal.values[psf.object$sink.nodes,]) ),
             beside=TRUE, col=bar.colors, names.arg=rep("",ncol(signal.values)),
             ylim=ylim, border=if (ncol(signal.values) < 80) "black" else NA )
    mtext(bquote("<log"[10] ~ "sinks>"), side=2, line=2.5, cex=1.5)
    
    ### Profile of mean node signals ###
    par(mar=c(6,7,0,5))
    
    bar.coords <- barplot( log10( colMeans(signal.values) ),
                           beside=TRUE, names.arg=rep("",ncol(signal.values)),
                           col=bar.colors, ylim=ylim, border=if (ncol(signal.values) < 80) "black" else NA )
    mtext(bquote("<log"[10] ~ "nodes>"), side=2, line=2.5, cex=1.5)
    
    if (ncol(signal.values)<100)
      text(bar.coords, par('usr')[3], labels=colnames(signal.values), srt=45, adj=c(1.1,1.1), xpd=TRUE)
  }
}


plot.psf.pathway <- function( pathway.object, signal.values=NULL, signal.values.lim, main="", highlight.genes=NULL, color.palette=NULL )
{
  width <- ncol(pathway.object$image)
  height <- nrow(pathway.object$image)
  max.xy <- max(width,height)
  
  par(mar = c(0, 0, 0, 0), mfrow=c(1,1) )
  plot(c(0-(max.xy-width)/2, width+(max.xy-width)/2), c(0-(max.xy-height)/2, height+(max.xy-height)/2), type = "n", xlab = "", ylab = "", xaxs = "i", yaxs = "i", axes=F)
  
  rasterImage(pathway.object$image, 0, 0, width, height, interpolate = F)
  
  title.node <- names( which( sapply( pathway.object$node.info, function(x) grepl("TITLE:", x$label ) ) ) )
  
  if(length(title.node)>0)
  {
    title.node <- pathway.object$node.info[[title.node]]
    
    rect( title.node$x-title.node$width*0.6, height-title.node$y+title.node$height*0.6,
          title.node$x+title.node$width*0.6, height-title.node$y-title.node$height*0.6,
          col="gray90", lwd=0.5 )
    
    text( title.node$x, height - title.node$y, main,
          cex = .5, col = "black", family="sans" )
  }
  
	
  gene.nodes <- names( which( sapply(pathway.object$node.info,"[[", "type" ) == "gene" ) )
  for( node in pathway.object$node.info[ gene.nodes ] )
  {
    node.col = "gray90"
    if( node$id %in% highlight.genes )
    {
      node.col = "olivedrab1"
      
    } else if( node$id %in% names(signal.values) )
    {
      node.col = color.palette(1000)[999*(signal.values[node$id]-signal.values.lim[1])/(signal.values.lim[2]-signal.values.lim[1])+1]
    }
    
    rect( node$x-node$width*0.5, height-node$y+node$height*0.5,
          node$x+node$width*0.5, height-node$y-node$height*0.5,
          col=node.col, lwd=0.5 )
    
    text( node$x, height - node$y, node$label,
          cex = ifelse( nchar(node$label)<=6,.45,.4), col = "black", family="sans" )
  }
  
  group.nodes <- names( which( sapply(pathway.object$node.info,"[[", "type" ) == "group" ) )
  for( node in pathway.object$node.info[ group.nodes ] )
  {
    node.col = "gray90"
    if( node$id %in% highlight.genes )
    {
      node.col = "olivedrab1"
      
    } else if( node$id %in% names(signal.values) )
    {
      node.col = color.palette(1000)[999*(signal.values[node$id]-signal.values.lim[1])/(signal.values.lim[2]-signal.values.lim[1])+1]
    }
    
    for( comp in node$components )
    {
      if(comp$type=="gene")
      {
        rect( comp$x-comp$width*0.5, height-comp$y+comp$height*0.5,
              comp$x+comp$width*0.5, height-comp$y-comp$height*0.5,
              col=node.col, lwd=0.5 )
        
        text( comp$x, height - comp$y, comp$label,
              cex = ifelse( nchar(comp$label)<=6,.45,.4), col = "black", family="sans" )
        
      } else if( comp$type=="compound")
      {
        circle( comp$x, height-comp$y, comp$width*0.75,
                col=node.col, lwd=0.5 )
      }
    }
  }
  
  compound.nodes <- names( which( sapply(pathway.object$node.info,"[[", "type" ) == "compound" ) )
  for( node in pathway.object$node.info[ compound.nodes ] )
    if( node$id %in% names(signal.values) )
    {
      node.col = color.palette(1000)[999*(signal.values[node$id]-signal.values.lim[1])/(signal.values.lim[2]-signal.values.lim[1])+1]

      circle( node$x, height-node$y, node$width*0.75,
              col=node.col, lwd=0.5 )
    }
  
  map.nodes <- names( which( sapply(pathway.object$node.info,"[[", "type" ) == "map" ) )
  map.nodes <- setdiff( map.nodes, title.node )
  
  for( node in pathway.object$node.info[ map.nodes ] )
  {
    rect( node$x-node$width*0.5-1, height-node$y+node$height*0.5+1,
          node$x+node$width*0.5+2, height-node$y-node$height*0.5-3,
          col = "gray95", lwd=0.5 )
    
    text( node$x, height - node$y, node$label,
          cex = 1-0.02*nchar(node$label), col = "black", family="sans" )
  }
  
}



psf.report.sheets <- function(env, psf.results, output.path, bar.colors )
{
  if( !setequal( rownames(psf.results[[1]]), names(kegg.collection[[1]]$node.info ) ) )
  {
    util.warn("PSF results outdated, so reports are skipped. Update oposSOM and re-run PSF calculation.")
    return()
  }
  
  util.info("Writing: PSF report sheets")
  progressbar <- newProgressBar(min = 0, max = length(kegg.collection)); cat("\r")
  
  for(pw in names(psf.results) )
  {
    filename <- file.path(output.path, paste( make.names(pw),".pdf",sep="") )
    pdf(filename, 29.7/2.54, 21/2.54)

    plot.psf.titlepage( env, psf.object=kegg.collection[[pw]], signal.values=psf.results[[pw]], bar.colors )
                        

    par(mfrow=c(1,1),mar=c(0,0,0,0))

    node.genes <- lapply(kegg.collection[[pw]]$node.info,"[[","ens.ids")
    node.genes <- lapply(node.genes,function(x) unique(env$gene.info$ensembl.mapping[which(env$gene.info$ensembl.mapping$ensembl_gene_id %in% x),1] ) )
    node.genes <- names( which( sapply(node.genes,function(x) any(!is.na(x)) ) ) )
 
    plot.psf.pathway( pathway.object=kegg.collection[[pw]], highlight.genes=node.genes, main=paste(pw,"\ngenes with data" ) )
 

    par(mfrow=c(1,1),mar=c(0,0,0,0))
    
    sink.genes <- kegg.collection[[pw]]$sink.nodes
    
    plot.psf.pathway( pathway.object=kegg.collection[[pw]], highlight.genes=sink.genes, main=paste(pw,"\nsink nodes") )
 
 
    for( m in 1:ncol(psf.results[[pw]]) )
    {
      signal.values = log10( psf.results[[pw]][,m] )

      plot.psf.pathway( pathway.object=kegg.collection[[pw]], signal.values = signal.values,
                        signal.values.lim = c(-1,1)*max( abs( log10( psf.results[[pw]] ) ) ),
                        main=paste(pw,"\n",colnames(psf.results[[pw]])[m] ),
                        color.palette=env$color.palette.heatmaps)
    }

    dev.off()

    setTxtProgressBar( progressbar, progressbar$getVal()+1 )
  }
  progressbar$kill()    

}

