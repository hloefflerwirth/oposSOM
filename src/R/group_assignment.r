pipeline.groupAssignment <- function()
{
  ### auto-assign PAT groups
 
  proto.pat.classification <- function(spotdata,proto.pat.data)
  {
    proto.pat.dist <- apply( proto.pat.data, 2, function(x) sum((x-spotdata)^2) )
    return( names(which.min(proto.pat.dist)) )
  }  
  
  if( length(unique(group.labels)==1) && unique(group.labels)[1] == "auto" )
  {
    util.info("Auto-assign sample groups (PAT-groups)")
    
    spot.list <- get(paste("spot.list.",preferences$standard.spot.modules,sep=""))
    
    pat.labels.sorted <- names( sort(table(pat.labels),decreasing=TRUE) )
    pat.labels.sorted <- pat.labels.sorted[which(pat.labels.sorted!="none")]
    
    withinss <- c()
    for( k in 2:min(100,length(pat.labels.sorted)) )
    {
      groups <- rep("",length(pat.labels))
      names(groups) <- names(pat.labels)
      
      groups[ which( pat.labels %in% pat.labels.sorted[1:k] ) ] <-
        pat.labels[ which( pat.labels %in% pat.labels.sorted[1:k] ) ]
      
      proto.pat.labels <- groups[ which(groups!="") ]
      proto.pat.data <- do.call( cbind, by( t(spot.list$spotdata[,names(proto.pat.labels)]), proto.pat.labels, colMeans ) )
      
      for( i in names(which(groups=="")) )
      {
        groups[i] <- proto.pat.classification(spot.list$spotdata[,i],proto.pat.data)
      }
      
      centroids <- sapply( unique(groups), function(x) rowMeans( spot.list$spotdata[,which(groups==x),drop=FALSE] )    )  
      withinss[as.character(k)] <- sum( sapply( 1:length(groups), function(i) sum((spot.list$spotdata[,i] - centroids[,groups[i]])^2) ) )
    }  
   
    
    x <- as.numeric(names(withinss))
    y <- withinss
    
    m <- (y[length(y)]-y[1]) / (x[length(x)]-x[1])
    n.ref <- y[1] - m * x[1]  

    # find optimal cluster number    
    opt.cluster.number <- length(withinss)
    for( k in 2:length(withinss) )
    {
      n.k <- y[as.character(k)] - m * x[k-1]
      if( all( n.k+m*x <= withinss ) )
      {
        opt.cluster.number <- k
        break        
      }
    } 
    n.opt <- y[as.character(opt.cluster.number)] - m * x[opt.cluster.number-1]
    
    
    
    opt.cluster.number <- max( min( opt.cluster.number + preferences$adjust.autogroup.number, ncol(indata) ), 1 )
    
    group.labels <<- rep("",length(pat.labels))
    names(group.labels) <<- names(pat.labels)
    
    group.labels[ which( pat.labels %in% pat.labels.sorted[1:opt.cluster.number] ) ] <<-
      pat.labels[ which( pat.labels %in% pat.labels.sorted[1:opt.cluster.number] ) ]
    
    proto.pat.labels <- group.labels[ which(group.labels!="") ]
    proto.pat.data <- do.call( cbind, by( t(spot.list$spotdata[,names(proto.pat.labels)]), proto.pat.labels, colMeans ) )
    
    for( i in names(which(group.labels=="")) )
    {
      group.labels[i] <<- proto.pat.classification(spot.list$spotdata[,i],proto.pat.data)
    }
    
    
    if(preferences$activated.modules$reporting)
    {  
      dir.create(paste(files.name, "- Results/Summary Sheets - Groups"), showWarnings=FALSE)
      
      filename <- file.path(paste(files.name, "- Results"),
                            "Summary Sheets - Groups",
                            "PAT-groups assignment.pdf")
      
      util.info("Writing:", filename)
      pdf(filename, 29.7/2.54, 21/2.54)
      
      x.coords <- barplot( withinss, col="gray90",ylim=range(y)*c(0.9,1.1), main="SSE",xlab="k",ylab="SSE", xpd=FALSE )
        box()
        lines( c(x.coords[1],x.coords[length(x)]), n.opt+m*c(x[1],x[length(x)]), col="red", xpd=FALSE )
        lines( c(x.coords[1],x.coords[length(x)]), n.ref+m*c(x[1],x[length(x)]), col="red", xpd=FALSE )
        abline(v=x.coords[opt.cluster.number-1],col="blue3",xpd=FALSE)
        abline(h=withinss[opt.cluster.number-1],col="blue3",xpd=FALSE)
      
      barplot( table(group.labels), main="PAT group frequency" )
      
      dev.off()
    }
    
    
    # sort data objects

    o <- order(group.labels)    
   
    pat.labels <<- pat.labels[o]
    group.labels <<- group.labels[o]
    group.colors <<- rep("", ncol(indata))
    for (i in seq_along(unique(group.labels)))
    {
      group.colors[which(group.labels == unique(group.labels)[i])] <<-
        colorRampPalette(c("blue3", "blue", "lightblue", "green2", "gold", "red", "red3"))(length(unique(group.labels)))[i]
    }
    names(group.colors) <<- names(group.labels)
    groupwise.group.colors <<- group.colors[match(unique(group.labels), group.labels)]
    names(groupwise.group.colors) <<- unique(group.labels)
    
    t.g.m <<- t.g.m[,o]
    p.g.m <<- p.g.m[,o]
    fdr.g.m <<- fdr.g.m[,o]
    Fdr.g.m <<- Fdr.g.m[,o]
    WAD.g.m <<- WAD.g.m[,o]
    n.0.m <<- n.0.m[o]
    perc.DE.m <<- perc.DE.m[o]    
    t.m <<- t.m[,o]
    p.m <<- p.m[,o]
    
    indata <<- indata[,o]
    metadata <<- metadata[,o]
    indata.sample.mean <<- indata.sample.mean[o]

    spot.list.samples <<- spot.list.samples[o]
    util.call(pipeline.detectSpotsIntegral,environment())
  }
  
  
  ### calculate group silouette
  
  if (length(unique(group.labels)) < 2)
  {
    group.silhouette.coef <<- rep(100, ncol(indata))
    return()
  }


  PCM <- cor( metadata )
  diag(PCM) <- NA
  
  group.silhouette.coef <<- sapply( seq(ncol(metadata)), function(i)
  {
    mean.group.correlations <- tapply( PCM[,i], group.labels, mean, na.rm=TRUE )
    
    cor.m.A <- mean.group.correlations[ group.labels[i] ]
    cor.m.B <- max( mean.group.correlations[ which( names(mean.group.correlations) != group.labels[i] ) ] )
    
    return( cor.m.A - cor.m.B )
  } )
  names(group.silhouette.coef) <<- colnames(metadata)

  group.silhouette.coef[which(is.nan(group.silhouette.coef))] <<- 0

}
