pipeline.summarySheetsPATs <- function()
{
  util.info("Processing PAT-centered Analyses")
  
  dir.create(paste(files.name, "- Results/Summary Sheets - PATs"), showWarnings=FALSE)
  
  pat.metadata <- do.call(cbind, by(t(metadata), pat.labels, colMeans))


  filename <- file.path(paste(files.name, "- Results"),
                        "Summary Sheets - PATs",
                        "PAT report.pdf")

  util.info("Writing:", filename)
  pdf(filename, 29.7/2.54, 21/2.54, useDingbats=FALSE)

  
  ### PAT clustering

  if( ncol(pat.metadata) > 2 )
  {
    hc <- hclust(dist(t(pat.metadata)))
  
    par(mar=c(8,6,4,5))
    plot(hc, main="PAT clustering", xlab="", sub="", hang=0.05 )
    
    cex.sample.portraits <- 0.025
    
    par(new=TRUE,mar=c(1,par("mar")[-1]))
    frame()
  
    for (i in seq_along(hc$order))
    {
      pat <- colnames(pat.metadata)[ hc$order[i] ]
                      
      m <- matrix(pat.metadata[,pat],
                  preferences$dim.1stLvlSom, preferences$dim.1stLvlSom)
  
      if (max(m) - min(m) != 0)  m <- 1 + (m - min(m)) / (max(m) - min(m)) * 999
  
      m <- cbind(apply(m, 1, function(x){x}))[nrow(m):1,]
      pix <- pixmapIndexed(m , col = color.palette.portraits(1000), cellres=10)
  
      addlogo(pix, (i-1)/(ncol(pat.metadata)-1)+cex.sample.portraits*c(-1,1), 0.05+0.08*(i%%2)+cex.sample.portraits*c(-1.2,1.2))
    }
  }

  ### PAT portraits
  
  layout(matrix(1:16,ncol=4,byrow=TRUE))
  
  for( pat in names(sort(table(pat.labels),decreasing=TRUE)) )
  {
       
    par(mar=c(0,0,0,0))
    frame()
      text(0.24, 0.94, pat, cex=3, adj=0)
  
    par(new=TRUE,mar=c(2,6,3,2))
    image(matrix(pat.metadata[,pat], preferences$dim.1stLvlSom, preferences$dim.1stLvlSom),
        axes=FALSE, col=color.palette.portraits(1000) )
        box()
         
    par(mar=c(0.5,0,4.5,2))  
    pat.samples <- names(which(pat.labels==pat))
    pie.table <- table(group.labels[pat.samples])[unique(group.labels)]
    pie.table <- pie.table[which(!is.na(pie.table))]
    pie(pie.table,col=groupwise.group.colors[names(pie.table)])

    mtext(paste("# samples:",sum(pat.labels==pat)),3)       
  }

  
  if( length(unique(pat.labels)) > 1 && length(unique(group.labels)) > 1 )
  {  
    ### Group to PAT relation barplots
    ##
    
    pat.group.assoc <- tapply(group.labels, pat.labels, function(x)
    {
      ret <- table(x)[unique(group.labels)]
      names(ret) <- unique(group.labels)
      return(ret)
    })
    pat.group.assoc <- do.call(cbind,pat.group.assoc)
    pat.group.assoc[which(is.na(pat.group.assoc))] <- 0
    
    layout(matrix(c(1, 2), 1, 2), c(4, 1), 1)
    par(mar=c(5, 5, 4, 2))
    barplot( pat.group.assoc[,ncol(pat.group.assoc):1],
            main="Association of PATs and groups", cex.main=1.5, cex.axis=1.4, cex.names=1,
            col=groupwise.group.colors, horiz=TRUE, las=1)
    
    par(mar=c(5, 1, 4, 2))
    plot(0, main="", cex.main=1, type="n", axes=FALSE, xlab="",
         ylab="", xlim=c(0,1), ylim=c(0,1), xaxs="i", yaxs="i")
  
    legend("topleft",legend=unique(group.labels),
         col=groupwise.group.colors, pch=15, pt.cex=1.4, title="groups")
  
    ##
    
    pat.group.assoc <- tapply(pat.labels, group.labels, function(x)
    {
      ret <- table(x)[sort(unique(pat.labels))]
      names(ret) <- sort(unique(pat.labels))
      return(ret)
    })[unique(group.labels)] 
    pat.group.assoc <- do.call(cbind,pat.group.assoc)
    pat.group.assoc[which(is.na(pat.group.assoc))] <- 0
  
    par(mfrow=c(1,1), mar=c(5, 8, 4, 2)) 
    y.coords <- barplot( pat.group.assoc[,ncol(pat.group.assoc):1],
             main="Association of PATs and groups", cex.main=1.5, cex.axis=1.4, cex.names=1,
             horiz=TRUE, las=1) 
    for( i in seq_along(unique(group.labels)) )
    {
      print.labels <- which( pat.group.assoc[,i] > max(colSums(pat.group.assoc,na.rm=TRUE))*0.05 )
      for( ii in seq_along(print.labels) )
      {
        text( x=sum( pat.group.assoc[1:(print.labels[ii]),i] ) - pat.group.assoc[print.labels[ii],i]/2,
              y=rev(y.coords)[i],
              names(print.labels)[ii] )
      }
    }
    
  }
  
  dev.off()
}
