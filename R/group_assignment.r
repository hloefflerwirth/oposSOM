pipeline.groupAssignment <- function()
{
  ### auto-assign PAT groups
 
  if( length(unique(group.labels)==1) && unique(group.labels)[1] == "auto" && any(pat.labels!="none")  )
  {
    util.info("Auto-assign sample groups (PAT-groups)")
    
    spot.list <- get(paste("spot.list.",preferences$standard.spot.modules,sep=""))
    pat.spotdata <- do.call( cbind, by( t(spot.list$spotdata), pat.labels, colMeans ) )
    sample.pat.distances <- apply( spot.list$spotdata, 2, function(x) apply( pat.spotdata, 2, function(y) sum((x-y)^2) ) )
    
    pat.labels.sorted <- names( sort(table(pat.labels),decreasing=TRUE) )
    pat.labels.sorted <- pat.labels.sorted[which(pat.labels.sorted!="none")[1]]
    pat.labels.test <- unique(pat.labels)[ which( !unique(pat.labels) %in% c( pat.labels.sorted, "none" ) ) ]
    
    if( length(pat.labels.test) <= 1 )
    {
      util.warn("Too few different PATs for clustering. Using all PATs as classes.")
      
      # assign new groups 
      
      group.labels <<- apply( sample.pat.distances[pat.labels.sorted[1:length(pat.labels.sorted)],,drop=FALSE], 2, function(x) names(x)[which.min(x)] )
      group.labels <<- paste( group.labels, "*" )
      names(group.labels) <<- colnames(indata)
      
    } else
    {
      withinss <- c()
      for( k in 2:min(100,length(unique(pat.labels))-1) )
      {
        test.pat.sse <- sapply( pat.labels.test, function(test.pat)
        {
          groups <- apply( sample.pat.distances[c(pat.labels.sorted,test.pat),], 2, function(x) names(x)[which.min(x)] )
          centroids <- pat.spotdata[,c(pat.labels.sorted,test.pat)]
          
          return( sum( sapply( names(groups), function(i) sum((spot.list$spotdata[,i] - centroids[,groups[i]])^2) ) ) )
        })
  
        pat.labels.sorted <- c(pat.labels.sorted, names(test.pat.sse)[which.min(test.pat.sse)] )
        pat.labels.test <- pat.labels.test[-which(pat.labels.test==names(test.pat.sse)[which.min(test.pat.sse)])]
        
        withinss[as.character(k)] <- min(test.pat.sse)
      }  
      
      x <- as.numeric(names(withinss))
      y <- withinss
      
      m <- if(length(y)>1) (y[length(y)]-y[1]) / (x[length(x)]-x[1]) else 0
      n.ref <- y[1] - m * x[1] 
      
      # find optimal cluster number 
      
      k.higher.sse <- sapply( 2:length(y), function(k) 
      {
        n.k <- y[as.character(k)] - m * x[k-1]
        sum( n.k+m*x <= y )
      })
      
      opt.cluster.number <- max( min( (which.max(k.higher.sse)+1) + preferences$adjust.autogroup.number, ncol(indata) ), 1 )
      n.opt <- y[as.character(which.max(k.higher.sse)+1)] - m * x[which.max(k.higher.sse)]
  
      # assign new groups
      
      group.labels <<- apply( sample.pat.distances[pat.labels.sorted[1:opt.cluster.number],,drop=FALSE], 2, function(x) names(x)[which.min(x)] )
      group.labels <<- paste( group.labels, "*" )
      names(group.labels) <<- colnames(indata)
      
      
      # plot SSE curve
      
      if(preferences$activated.modules$reporting)
      {  
        dir.create(paste(files.name, "- Results/Summary Sheets - Groups"), showWarnings=FALSE)
        
        filename <- file.path(paste(files.name, "- Results"),
                              "Summary Sheets - Groups",
                              "PAT-groups assignment.pdf")
        
        util.info("Writing:", filename)
        pdf(filename, 29.7/2.54, 21/2.54, useDingbats=FALSE)
        
        x.coords <- barplot( withinss, col="gray90",ylim=range(y)*c(0.9,1.1), main="SSE",xlab="k",ylab="SSE", xpd=FALSE )
          box()
          lines( c(x.coords[1],x.coords[length(x)]), n.opt+m*c(x[1],x[length(x)]), col="red", xpd=FALSE )
          lines( c(x.coords[1],x.coords[length(x)]), n.ref+m*c(x[1],x[length(x)]), col="red", xpd=FALSE )
          abline(v=x.coords[opt.cluster.number-1],col="blue3",xpd=FALSE)
          abline(h=y[opt.cluster.number-1],col="blue3",xpd=FALSE)
        
        barplot( table(group.labels), main="PAT group frequency" )
        
        dev.off()
      }
    
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
