pipeline.genesetProfilesAndMaps <- function(env)
{
  dirname <- file.path(paste(env$files.name, "- Results"), "Geneset Analysis")
  util.info("Writing:", file.path(dirname, "*.{csv,pdf}"))
	progressbar <-newProgressBar(min = 0, max = nrow(env$samples.GSZ.scores)); cat("\r")
  
  ylim <- quantile(env$samples.GSZ.scores,c(0.01,0.99))
  off.thres <- -sd(env$samples.GSZ.scores)
  on.thres <- sd(env$samples.GSZ.scores)
  
  
  for (i in 1:nrow(env$samples.GSZ.scores))
  {
    filename.prefix <- substring(make.names(names(env$gs.def.list)[i]), 1, 80)
    pdf(file.path(dirname, paste(filename.prefix, ".pdf",sep="")), 29.7/2.54, 21/2.54, useDingbats=FALSE)
     
    #### Geneset Profile + Heatmap
    layout(matrix(c(1,2,3),ncol=1,byrow=TRUE),heights=c(1.5,0.5,4))
    
    set.genes <- unique( env$gene.info$ensembl.mapping[which(env$gene.info$ensembl.mapping$ensembl_gene_id %in% env$gs.def.list[[i]]$Genes),1] )
    
    gs.indata <- env$indata[set.genes,,drop=F]
    sig.genes <- which( apply(gs.indata,1,sd)>sd(gs.indata) )
    if(length(sig.genes)==0) sig.genes <- order(apply(gs.indata,1,sd),decreasing=TRUE)[1]

    if( ncol(env$indata) * length(sig.genes) < 20000 )
    {
      gs.indata <- gs.indata[ sig.genes, ,drop=FALSE]
      rownames(gs.indata) <- env$gene.info$names[rownames(gs.indata)]
      
      o.genes <- if(length(sig.genes)>1) hclust(dist(gs.indata))$order else 1
      o.samples <- order(env$samples.GSZ.scores[i,])
      
      
      offset <- min(env$samples.GSZ.scores[i,])
      
      par(mar=c(0,10,5,16))
      barplot( env$samples.GSZ.scores[i,o.samples]-offset, ylab="", xlab="", ylim=range(env$samples.GSZ.scores[i,])-offset, xaxt="n", xaxs="i",
                col=env$group.colors[o.samples], xpd=TRUE, space=c(0,0), offset=0, axes=FALSE, border=if (ncol(env$indata) < 100) "black" else NA )
        abline( h=0-offset, lty=2, col="gray" )
        axis( 2, at=seq(from=0,to=(max(env$samples.GSZ.scores[i,])-offset),length.out=4), labels=round(seq(from=offset,to=max(env$samples.GSZ.scores[i,]),length.out=4),2),las=2, cex.axis=1.4)
        mtext( names(env$gs.def.list[i]), side=3, cex=1.5, line=1 )
        mtext("GSZ", side=2, line=4.5, cex=1.25)
        
  
      par(new=TRUE,mar=c(0,0,0,0))
      frame()
        legend(x=0.86,y=0.7,names(env$groupwise.group.colors),text.col=env$groupwise.group.colors)
  
      par(mar=c(0.5,10,0.5,16))
      image( matrix(1:ncol(env$indata),ncol(env$indata),1), col=env$group.colors[o.samples], axes=FALSE )
        box()
      
      par(mar=c(10,10,0,16))
      image(x=1:ncol(gs.indata), y=1:nrow(gs.indata),z=t(gs.indata[o.genes,o.samples,drop=FALSE]),col=env$color.palette.heatmaps(1000), axes=FALSE,zlim=max(max(gs.indata),-min(gs.indata))*c(-1,1),xlab="", ylab="")
        box()
        axis(2, 1:nrow(gs.indata), labels=rownames(gs.indata)[o.genes], las=2, line=-0.5, tick=0, cex.axis=min( 1.5, max( 0.5, 2-nrow(gs.indata) /67 ) ) )
      
      par(new=TRUE,mar=c(40,74,0,1))
      image(matrix(c(1:1000), 1000, 1), axes=FALSE, col=env$color.palette.heatmaps(1000))
        box()
        axis( 1, at=c(0,1), labels=round(range(gs.indata),1), cex.axis=1.4 )
        mtext( bquote(Delta ~ "e"), side=1, cex=1.25, line=1 )  
    }
      
      
    #### Geneset Profiles + Mapping
    layout(matrix(c(1,2,3,0),ncol=2),widths=c(1,0.5), heights=c(1,0.08))
      
    # barplot
    par(mar=c(14,4,4,1))
    barplot(env$samples.GSZ.scores[i,], beside=TRUE,
            las=2, cex.names=1.2, cex.axis=1.4, col=env$group.colors, 
            ylim=ylim, border=if (ncol(env$indata) < 100) "black" else NA,
            names.arg=rep("",ncol(env$indata)))

      abline( h=c(off.thres,on.thres,0), lty=2)
      mtext( names(env$gs.def.list[i]), side=3, cex=1.5, line=1 )
      mtext("GSZ", side=2, line=2.5, cex=1.25)
    
    # barcode
    par(mar=c(1,5.5,0.1,2.7))
    
    col <- rep("gray60",ncol(env$indata))
    col[which(env$samples.GSZ.scores[i,]<off.thres)] <- "white"
    col[which(env$samples.GSZ.scores[i,]>on.thres)] <- "black"
    
    image( matrix(1:ncol(env$indata),ncol(env$indata),1), col=col, axes=FALSE )
     box()
    
    par(new=TRUE, mar=c(1,1,0.1,4))
    frame()
     legend("left",c("on","-","off"),fill=c("black","gray60","white"),border="black",bty="n")
    
    # population map
    par(mar=c(16, 0, 6, 3.5))
     
    n.map <- matrix(0,env$preferences$dim.1stLvlSom,env$preferences$dim.1stLvlSom)
    gs.nodes <- env$som.result$feature.BMU[set.genes]
    n.map[as.numeric(names(table(gs.nodes)))] <- table(gs.nodes)
    n.map[which(n.map==0)] <- NA
    n.map <- matrix(n.map, env$preferences$dim.1stLvlSom)
    
    lim <- c(1,env$preferences$dim.1stLvlSom) + env$preferences$dim.1stLvlSom * 0.01 * c(-1, 1)
    colr <- env$color.palette.heatmaps(1000)[(na.omit(as.vector(n.map)) - min(n.map,na.rm=TRUE)) /
                           max(1, (max(n.map,na.rm=TRUE) - min(n.map,na.rm=TRUE))) *
                           999 + 1]
    
    plot(which(!is.na(n.map), arr.ind=TRUE), xlim=lim, ylim=lim, pch=16, axes=FALSE,
        xlab="",ylab="", xaxs="i", yaxs="i", col=colr,
        cex=0.5 + na.omit(as.vector(n.map)) / max(n.map,na.rm=TRUE) * 2.8)
    
      title(sub=paste("# features =", length(set.genes)),line=0)
      box()
    
    par(new=TRUE, mar=c(32,21.5,6,0.5))
    image(matrix(1:100, 1, 100), col = env$color.palette.heatmaps(1000), axes=FALSE)
      axis(2, at=c(0,1), c(min(n.map,na.rm=TRUE), max(n.map,na.rm=TRUE)), las=2, tick=FALSE, pos=-0.36)
      box()

      
    #### Geneset Profiles - Group boxplots
    layout(1)
    par(mar=c(12,6,4,5))
    
    mean.boxes <- by(env$samples.GSZ.scores[i,], env$group.labels, c)[unique(env$group.labels)]
    boxplot(mean.boxes, col=env$groupwise.group.colors, ylim=ylim+sum(abs(ylim))*0.1*c(-1,1), axes=FALSE, yaxs="i")
      box()
      axis(1, seq_along(unique(env$group.labels)), unique(env$group.labels), las=2, tick=FALSE)
      axis(2, las=2)

      abline(h=0, lty=2)

      mtext( names(env$gs.def.list[i]), side=3, cex=1.5, line=1 )
      mtext("GSZ", side=2, line=2.5, cex=1.25)

    dev.off()
  

    ##### CSV sheets

    out <- data.frame(AffyID=names(env$gene.info$ids[set.genes]),
                      EnsemblID=env$gene.info$ids[set.genes],
                      Metagene=env$gene.info$coordinates[set.genes],
                      Max.expression.sample=colnames(env$indata)[apply(env$indata[set.genes, ,drop=FALSE], 1, which.max)],
                      GeneSymbol=env$gene.info$names[set.genes],
                      Description=env$gene.info$descriptions[set.genes])

    env$csv.function(out, file.path(dirname, paste(filename.prefix, ".csv", sep="")), row.names=FALSE)

    setTxtProgressBar( progressbar, progressbar$getVal()+1 )
  }

  progressbar$kill()
}
