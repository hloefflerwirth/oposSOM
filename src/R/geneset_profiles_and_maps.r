pipeline.genesetProfilesAndMaps <- function()
{
  progress.current <- 0
  progress.max <- nrow(samples.GSZ.scores) + length(gs.def.list)
  util.progress(progress.current, progress.max)


#   if (preferences$geneset.analysis.exact)
#   {
#     fdr.threshold <- 0.1
# 
#     GSZs <-
#       cbind(as.vector(unlist(sapply(spot.list.samples, function(x) { x$GSZ.score[names(gs.def.list)] }))),
#             as.vector(unlist(sapply(spot.list.samples, function(x) { x$GSZ.p.value[names(gs.def.list)] }))))
# 
#     GSZs <- GSZs[which(!is.na(GSZs[,2])),]
# 
#     fdrtool.result <- fdrtool(GSZs[,2], statistic="pvalue", verbose=FALSE, plot=FALSE)
#     fdr.significant.spot.list.samples <- length(which(fdrtool.result$lfdr < fdr.threshold))
#     fdr.gsz.threshold <- sort(GSZs[,1], decreasing=TRUE)[fdr.significant.spot.list.samples]
#   } else
#   {
#     fdr.gsz.threshold <- 0
#   }


  ## Geneset Profiles over Samples
  dirname <- file.path(paste(files.name, "- Results"), "Geneset Analysis")
  util.info("Writing:", file.path(dirname, "*.{csv,pdf}"))

  for (i in 1:nrow(samples.GSZ.scores))
  {
    filename.prefix <- substring(make.names(names(gs.def.list)[i]), 1, 100)
    pdf(file.path(dirname, paste(filename.prefix, "profile.pdf")), 29.7/2.54, 21/2.54)
     
    
    dens <- density(samples.GSZ.scores[i,],adjust=1)
    
    dens.dev <- dens$y[1:511] - dens$y[2:512]
    dens.extrema <- which( sign(dens.dev[1:511]) != sign(dens.dev[2:512]) )
    dens.off.peak <- min(dens.extrema)
    dens.on.peak <- max(dens.extrema)
    
    if( dens.on.peak != dens.off.peak )
    {    
      dens.off.thres <- dens$x[ min( 2*dens.off.peak, 512 ) ]
      dens.on.thres <- dens$x[ max( dens.on.peak-(512-dens.on.peak), 1 ) ]
            
      if( sum(samples.GSZ.scores[i,]<dens.off.thres) > 3 && sum(samples.GSZ.scores[i,]>dens.on.thres) > 3 )
      {
        dens.left <- dens
        dens.left$y[ dens.off.peak:min( (2*dens.off.peak-1), 512 ) ] <- rev( dens.left$y[ 1:length(dens.off.peak:min( (2*dens.off.peak-1), 512 )) ] )
        dens.left$y[ min( 2*dens.off.peak, 512 ):512 ] <- 0
        dens.right <- dens
        dens.right$y[ max( (dens.on.peak-(512-dens.on.peak)), 1 ):dens.on.peak ] <- rev( dens.right$y[ dens.on.peak:( dens.on.peak + length(max( (dens.on.peak-(512-dens.on.peak)), 1 ):dens.on.peak) - 1 ) ] )
        dens.right$y[ 1:max( (dens.on.peak-(512-dens.on.peak)-1), 1 ) ] <- 0    
        
        dens.off.sd <- sd(samples.GSZ.scores[i,][which(samples.GSZ.scores[i,]<dens.off.thres)])
        dens.on.sd <- sd(samples.GSZ.scores[i,][which(samples.GSZ.scores[i,]>dens.on.thres)])      
        
        dens.off.thres <- dens$x[dens.off.peak]+dens.off.sd
        dens.on.thres <- dens$x[dens.on.peak]-dens.on.sd
                
        if ( dens.on.thres < dens.off.thres )
        {
          h <- dens.on.thres
          dens.on.thres <- dens.off.thres
          dens.off.thres <- h
        }
        
      } else
      {      
        dens.on.peak <- dens.off.peak
        dens.off.thres <- -sd(samples.GSZ.scores[i,])
        dens.on.thres <- sd(samples.GSZ.scores[i,])
      }
            
    } else
    {      
      dens.off.thres <- -sd(samples.GSZ.scores[i,])
      dens.on.thres <- sd(samples.GSZ.scores[i,])
    }
    
    
    
    layout(matrix(c(1,3,2,2),ncol=2),widths=c(1,0.2), heights=c(1,0.08))
       
    
    # barplot
    ylim <- c(-10, 20)
    
    par(mar=c(14,6,4,2))
    bar.coords <- barplot(samples.GSZ.scores[i,], beside=TRUE,
                          main=rownames(samples.GSZ.scores)[i], las=2,
                          cex.names=1.2, col=group.colors, cex.main=1,
                          ylim=ylim,
                          border=if (ncol(indata) < 80) "black" else NA,
                          names.arg=if (ncol(indata)<100) colnames(indata) else rep("",ncol(indata)))
    
    abline(h=c(dens.off.thres,dens.on.thres), lty=2)
    
    
    mtext("GSZ", side=2, line=2.5, cex=1.5)
    
    
    
    # density
    ylim <- c(-21.8, 21.5)
    
    par(mar=c(5,0.1,1,2))
    plot(rev(dens$y), rev(dens$x), type="l", xlab="", ylab="", axes=FALSE, ylim=ylim, lwd=2, col="gray" )
    if( dens.on.peak != dens.off.peak )
    {  
      lines( rev(dens.left$y), rev(dens.left$x) )
      lines( rev(dens.right$y), rev(dens.right$x) )
    
      abline( h=c(dens$x[dens.off.peak],dens$x[dens.on.peak]) )
      
      abline( h=dens$x[dens.off.peak]-dens.off.sd, lty=2, col="gray" )
      abline( h=dens$x[dens.on.peak]+dens.on.sd, lty=2, col="gray" )
    }
    
    abline(h=c(dens.off.thres,dens.on.thres),lty=2)
    
    box()
    
    
    
    # barcode
    par(mar=c(1,8,0.1,4))
    
    col = rep("gray60",ncol(indata))
    col[which(samples.GSZ.scores[i,]<dens.off.thres)] = "white"
    col[which(samples.GSZ.scores[i,]>dens.on.thres)] = "black"
    
    image( matrix(1:ncol(indata),ncol(indata),1), col=col, axes=FALSE )
    box()
    
    par(new=TRUE, mar=c(1,3,0.1,4))
    frame()
    legend("left",c("on","-","off"),fill=c("black","gray60","white"),border="black",bty="n")
    
    
    

    #################################################

    if( length(spot.list.overexpression$spots) > 2 )
    {
      par(mfrow=c(1,2))
      par(mar=c(5,6,4,5))
      
      spot.fisher.p <- -log10(sapply(spot.list.overexpression$spots, function(x)
      {
        x$Fisher.p[names(gs.def.list)[i]]
      }))
      names(spot.fisher.p) <- LETTERS[seq_along(spot.fisher.p)]
      spot.fisher.p <- pmin( spot.fisher.p, 16 )
  
      radarchart( as.data.frame( rbind( rep(10,length(spot.fisher.p) ), rep(0,length(spot.fisher.p) ), spot.fisher.p ) ),   						
                  pcol="gray50", pfcol="gray80",	title="Enrichment in overexpression spots", seg=2, pty=32, plwd = 3 )
      
      
      
      spot.fisher.p <- -log10(sapply(spot.list.kmeans$spots, function(x)
      {
        x$Fisher.p[names(gs.def.list)[i]]
      }))  
      names(spot.fisher.p) <- LETTERS[seq_along(spot.fisher.p)]
      spot.fisher.p <- pmin( spot.fisher.p, 16 )
  
      radarchart( as.data.frame( rbind( rep(10,length(spot.fisher.p) ), rep(0,length(spot.fisher.p) ), spot.fisher.p ) ),     					
                  pcol="gray50", pfcol="gray80",	title="Enrichment in kMeans clusters", seg=2, pty=32, plwd = 3 )
    }
      
    
    dev.off()

    #################################################

    genes <- names(gene.ids)[which(gene.ids %in% gs.def.list[[i]]$Genes)]

    out <- data.frame(AffyID=names(gene.ids[genes]),
                      EnsemblID=gene.ids[genes],
                      Metagene=gene.coordinates[genes],
                      Max.expression.sample=colnames(indata)[apply(indata[genes, ,drop=FALSE], 1, which.max)],
                      GeneSymbol=gene.names[genes],
                      Description=gene.descriptions[genes])

    write.csv2(out, file.path(dirname, paste(filename.prefix, ".csv", sep="")), row.names=FALSE)

    progress.current <- progress.current + 1
    util.progress(progress.current, progress.max)
  }

  ## Gensets Population Maps
  for (i in seq_along(gs.def.list))
  {
    filename <- paste(substring(make.names(names(gs.def.list)[i]), 1, 100), "map.pdf")
    pdf(file.path(dirname, filename), 21/2.54, 21/2.54)

    ### Population Map ###
    n.map <- matrix(0,preferences$dim.1stLvlSom,preferences$dim.1stLvlSom)
    gs.nodes <- som.nodes[names(gene.ids)[which(gene.ids %in% gs.def.list[[i]]$Genes)]]
    n.map[as.numeric(names(table(gs.nodes)))] <- table(gs.nodes)
    n.map[which(n.map==0)] <- NA
    n.map <- matrix(n.map, preferences$dim.1stLvlSom)

    x.test <-
      sum((apply(n.map, 1, sum, na.rm=TRUE) - sum(n.map, na.rm=TRUE)*(1/preferences$dim.1stLvlSom))^2 /
          sum(n.map, na.rm=TRUE) *
          (1/preferences$dim.1stLvlSom))

    y.test <-
      sum((apply(n.map, 2, sum, na.rm=TRUE) - sum(n.map, na.rm=TRUE)*(1/preferences$dim.1stLvlSom))^2 /
          sum(n.map, na.rm=TRUE) *
          (1/preferences$dim.1stLvlSom))

    chi.sq.p.value <- 1-pchisq(x.test+y.test, 1)

    par(mfrow=c(1, 1))
    par(mar=c(5, 6, 4, 5))

    lim <- c(1,preferences$dim.1stLvlSom) + preferences$dim.1stLvlSom * 0.01 * c(-1, 1)

    colr <- colramp(1000)[(na.omit(as.vector(n.map)) - min(n.map,na.rm=TRUE)) /
                          max(1, (max(n.map,na.rm=TRUE) - min(n.map,na.rm=TRUE))) *
                          999 + 1]

    plot(which(!is.na(n.map), arr.ind=TRUE), xlim=lim, ylim=lim, pch=16, axes=FALSE,
         xlab="",ylab="", xaxs="i", yaxs="i", main=names(gs.def.list)[i], col=colr,
         cex=0.5 + na.omit(as.vector(n.map)) / max(n.map,na.rm=TRUE) * 2.8)

    title(sub=paste("# features =", sum(gene.ids %in% gs.def.list[[i]]$Genes)),line=0)
    title(sub=paste("chi-square p = ",round(chi.sq.p.value,2)),line=1)
    box()

    par(new=TRUE, mar=c(1,0,0,0))
    layout(matrix(c(0, 0, 0, 0, 1, 0, 0, 0, 0), 3, 3), c(1, 0.05, 0.02), c(0.15, 0.3, 1))
    image(matrix(1:100, 1, 100), col = colramp(1000), axes=FALSE)
    axis(2, at=c(0,1), c(min(n.map,na.rm=TRUE), max(n.map,na.rm=TRUE)), las=2, tick=FALSE, pos=-0.5)
    box()

    ### Smooth Population Map ###
    if (any(!is.na(n.map)))
    {
      par(mfrow=c(1, 1))
      par(mar=c(5, 6, 4, 5))

      coords <-
        som.result$visual[which(rownames(indata) %in% names(gene.ids[which(gene.ids %in% gs.def.list[[i]]$Genes)])),
                          c(1,2), drop=FALSE]

      if (nrow(coords) == 1)
      {
        coords <- rbind(coords, coords)
      }

      suppressWarnings(smoothScatter(coords+1 ,
                                     main=names(gs.def.list)[i],
                                     cex.main=2,
                                     xlim=c(1,preferences$dim.1stLvlSom),
                                     ylim=c(1,preferences$dim.1stLvlSom),
                                     xaxs="i", yaxs="i", axes=FALSE, xlab="",
                                     ylab="", nrpoints=0))

      title(sub=paste("# features =", sum(gene.ids %in% gs.def.list[[i]]$Genes),
                      ", max =", max(n.map,na.rm=TRUE)),line=0)
    }

    dev.off()

    progress.current <- progress.current + 1
    util.progress(progress.current, progress.max)
  }

  util.progress.terminate()
}
