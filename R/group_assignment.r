pipeline.groupAssignment <- function(env)
{
  ### auto-assign PAT groups
 
  if( length(unique(env$group.labels)==1) && unique(env$group.labels)[1] == "auto" && any(env$pat.labels!="none")  )
  {
    util.info("Auto-assign sample groups (PAT-groups)")
    
    spot.list <- env[[paste("spot.list.", env$preferences$standard.spot.modules,sep="")]]
    pat.spotdata <- do.call( cbind, by( t(spot.list$spotdata), env$pat.labels, colMeans ) )
    sample.pat.distances <- apply( spot.list$spotdata, 2, function(x) apply( pat.spotdata, 2, function(y) sum((x-y)^2) ) )
    
    pat.labels.sorted <- names( sort(table(env$pat.labels),decreasing=TRUE) )
    pat.labels.sorted <- pat.labels.sorted[which(pat.labels.sorted!="none")[1]]
    pat.labels.test <- unique(env$pat.labels)[ which( !unique(env$pat.labels) %in% c( pat.labels.sorted, "none" ) ) ]
    
    if( length(pat.labels.test) <= 1 )
    {
      util.warn("Too few different PATs for clustering. Using all PATs as classes.")
      
      # assign new groups 
      
      env$group.labels <- apply( sample.pat.distances[pat.labels.sorted[1:length(pat.labels.sorted)],,drop=FALSE], 2, function(x) names(x)[which.min(x)] )
      env$group.labels <- paste( env$group.labels, "*" )
      names(env$group.labels) <- colnames(env$indata)
      
    } else
    {
      withinss <- c()
      for( k in 2:min(100,length(unique(env$pat.labels))-1) )
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
      
      opt.cluster.number <- max( min( (which.max(k.higher.sse)+1) + env$preferences$adjust.autogroup.number, ncol(env$indata) ), 1 )
      n.opt <- y[as.character(which.max(k.higher.sse)+1)] - m * x[which.max(k.higher.sse)]
  
      # assign new groups
      
      env$group.labels <- apply( sample.pat.distances[pat.labels.sorted[1:opt.cluster.number],,drop=FALSE], 2, function(x) names(x)[which.min(x)] )
      env$group.labels <- paste( env$group.labels, "*" )
      names(env$group.labels) <- colnames(env$indata)
      
      
      # plot SSE curve
      
      if(env$preferences$activated.modules$reporting)
      {  
        dir.create(paste(env$files.name, "- Results/Summary Sheets - Groups"), showWarnings=FALSE)
        
        filename <- file.path(paste(env$files.name, "- Results"),
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
        
        barplot( table(env$group.labels), main="PAT group frequency" )
        
        dev.off()
      }
    
    }  
      
    # sort data objects

    o <- order(env$group.labels)    
   
    env$pat.labels <- env$pat.labels[o]
    env$group.labels <- env$group.labels[o]
    env$group.colors <- rep("", ncol(env$indata))
    for (i in seq_along(unique(env$group.labels)))
    {
      env$group.colors[which(env$group.labels == unique(env$group.labels)[i])] <- color.palette.discrete(length(unique(env$group.labels)))[i]
    }
    names(env$group.colors) <- names(env$group.labels)
    env$groupwise.group.colors <- env$group.colors[match(unique(env$group.labels), env$group.labels)]
    names(env$groupwise.group.colors) <- unique(env$group.labels)
    
    env$p.g.m <- env$p.g.m[,o]
    env$fdr.g.m <- env$fdr.g.m[,o]
    env$n.0.m <- env$n.0.m[o]
    env$perc.DE.m <- env$perc.DE.m[o]    
    env$p.m <- env$p.m[,o]
    
    env$indata <- env$indata[,o]
    env$metadata <- env$metadata[,o]
    env$indata.sample.mean <- env$indata.sample.mean[o]

    env$spot.list.samples <- env$spot.list.samples[o]
    env <- pipeline.detectSpotsIntegral(env)
  }
  
  
  ### calculate group silouette
  
  if (length(unique(env$group.labels)) < 2)
  {
    env$group.silhouette.coef <- rep(100, ncol(env$indata))
    return()
  }


  PCM <- cor( env$metadata )
  diag(PCM) <- NA
  
  env$group.silhouette.coef <- sapply( seq(ncol(env$metadata)), function(i)
  {
    mean.group.correlations <- tapply( PCM[,i], env$group.labels, mean, na.rm=TRUE )
    
    cor.m.A <- mean.group.correlations[ env$group.labels[i] ]
    cor.m.B <- max( mean.group.correlations[ which( names(mean.group.correlations) != env$group.labels[i] ) ] )
    
    return( cor.m.A - cor.m.B )
  } )
  names(env$group.silhouette.coef) <- colnames(env$metadata)

  env$group.silhouette.coef[which(is.nan(env$group.silhouette.coef))] <- 0
  return(env)
}
