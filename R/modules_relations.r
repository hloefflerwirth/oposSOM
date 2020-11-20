modules.relations <- function(env, spot.list, main, path)
{
  sd.theshold <- sd(spot.list$spotdata)
  spotdata.binary <- spot.list$spotdata > sd.theshold
  mode(spotdata.binary) <- "numeric"

  pdf(path, 29.7/2.54, 21/2.54, useDingbats=FALSE)

  #### wTO network ####
  
  if ( length(spot.list$spots) > 2 )
  {
    # calculate wTO #
    adj_matrix <- cor(t(spot.list$spotdata),t(spot.list$spotdata))
    diag(adj_matrix) <- 0
    
    omega <- matrix(NA,
                    nrow=nrow(spot.list$spotdata),
                    ncol=nrow(spot.list$spotdata),
                    dimnames=list(names(spot.list$spots),names(spot.list$spots)))
    
    b <- adj_matrix %*% adj_matrix
    k_sum <- apply(adj_matrix, 2, function(x) { sum(abs(x)) })
    
    for (i in 1: nrow(spot.list$spotdata))
    {
      for (j in 1: nrow(spot.list$spotdata))
      {
        omega[i,j] <-
          (b[i,j] + adj_matrix[i,j]) /
          (min(k_sum[i],k_sum[j]) + 1 - abs(adj_matrix[i,j]))
      }
    }
    
    diag(omega) <- 0
    omega[which(abs(omega) < 0.25)] <- 0
    
    g <- graph.adjacency(omega, weighted=TRUE, mode="undirected")

    V(g)$label <- names(spot.list$spots)
    n.spot <- apply(spot.list$spotdata, 1, function(x) sum(x>sd.theshold) )
    
    if (max(n.spot) == min(n.spot))
    {
      V(g)$size <- 10
    } else
    {
      V(g)$size <- 8 * (n.spot-min(n.spot)) / (max(n.spot)-min(n.spot)) + 6
    }
    
    if (length(E(g)) > 0)
    {
      E(g)$color <- c("brown4","chocolate","rosybrown1","white","white","lightblue","cornflowerblue","dodgerblue4")[ cut( E(g)$weight, breaks=c(-Inf,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,Inf), labels=c(1:8) ) ]
      E(g)$width <- abs(E(g)$weight)
      
      if (max(E(g)$width) == min(E(g)$width))
      {
        E(g)$width <- 3
      } else
      {
        E(g)$width <- 5 * (E(g)$width-min(E(g)$width)) /
          (max(E(g)$width)-min(E(g)$width)) + 0.5
      }
    }
    
    # plot wTO network
    
    layout(matrix(c(1, 2), 1, 2), c(2, 1), 1)
    par(mar=c(5, 4, 4, 1.8))
    
    image(matrix(spot.list$overview.mask, env$preferences$dim.1stLvlSom), col="gray90", axes=FALSE, main="Module correlations (WTO algorithm)", cex.main=1.6)
      box()
    
    par(new=TRUE, mar=c(3.9,2.7,2.7,0.5))
    plot(g, layout=t(sapply(spot.list$spots, function(x) { x$position })),
        vertex.label.color="black" ,vertex.label.cex=1, vertex.color="grey",
        xlim=c(-1.0,1.0), ylim=c(-1.0,1.0))
  
    par(mar=c(5, 1, 4, 2))
    plot(0, type="n", axes=FALSE, xlab="", ylab="", xlim=c(0,1), ylim=c(0,1), xaxs="i", yaxs="i")
    legend("center", c("negative correlation",expression(paste("  ",omega," < -0.25")),expression(paste("  ",omega," < -0.50")),expression(paste("  ",omega," < -0.75")),
                          "positive correlation",expression(paste("  ",omega," > 0.25")),expression(paste("  ",omega," > 0.50")),expression(paste("  ",omega," > 0.75"))),
                          lty = c(1,1), lwd = 4,col=c("white","rosybrown1","chocolate","brown4","white","lightblue","cornflowerblue","dodgerblue4"), bty="n")
    box()
  }
 
  
  #### Baskets ####
  
  if ( sum( rowSums(spotdata.binary) > 0 ) <=1 )
  {
    util.warn("Skipped basket and group association analyses due to missing module activation.")
    return()
  }
  
  if ( length(spot.list$spots) > 2 )
  {
    module.group.rules <-
      do.call( rbind, lapply( rownames(spotdata.binary), function(mod.x)
      {
        samples.with.module <- names( which( spotdata.binary[mod.x,] == 1 ) )
        do.call( rbind, lapply( unique( env$group.labels[ samples.with.module ] ), function(gr.x)
        {
          count <- sum( env$group.labels[ samples.with.module ] == gr.x )
          data.frame( lhs = mod.x,
                      rhs = gr.x,
                      count = count,
                      support = count / ncol(env$indata),
                      confidence = count / length(samples.with.module),
                      stringsAsFactors = FALSE)
        }) )
      } ) )
    module.group.rules <- module.group.rules[ which( module.group.rules$support>0.01 &
                                                     module.group.rules$confidence>0.1 ), ]
    module.group.rules <- module.group.rules[ order( module.group.rules$confidence, decreasing=TRUE), ]
    module.group.rules <- module.group.rules[ !duplicated( module.group.rules$lhs ) ,] 
    
    
    module.module.rules <-
      do.call( rbind, lapply( rownames(spotdata.binary), function(mod.x.1)
      {
        do.call( rbind, lapply( setdiff( rownames(spotdata.binary), mod.x.1), function(mod.x.2)
        {
          count <- sum( colSums( spotdata.binary[c(mod.x.1,mod.x.2),] ) == 2 )
          data.frame( lhs = mod.x.1,
                      rhs = mod.x.2,
                      count = count,
                      support = count / ncol(env$indata),
                      confidence = count / sum( spotdata.binary[mod.x.1,] == 1 ),
                      stringsAsFactors = FALSE)
        }) )
      } ) )
    module.module.rules <- module.module.rules[ which( module.module.rules$support>0.01 &
                                                         module.module.rules$confidence>0.1 ), ]
    module.module.rules <- module.module.rules[ order( module.module.rules$confidence, decreasing=TRUE), ]
    module.module.rules <- module.module.rules[ which( !duplicated( apply(cbind(module.module.rules$lhs,module.module.rules$rhs),1,function(x) paste(sort(x),collapse=" ") ) ) ),]
    
  
    g <- graph.empty( length(spot.list$spots), directed = FALSE )
    V(g)$name <- names(spot.list$spots)
    g <- add_edges( g, as.vector( rbind( match( module.module.rules$lhs, V(g)$name ), match( module.module.rules$rhs, V(g)$name ) ) ) )
    V(g)$label <- names(spot.list$spots)
    V(g)$size <- 5
    V(g)$size[ match( module.group.rules$lhs, V(g)$name ) ] <- 15 * module.group.rules$confidence + 3
    V(g)$color <- "gray"
    V(g)$color[ match( module.group.rules$lhs, V(g)$name ) ] <- env$groupwise.group.colors[module.group.rules$rhs]
    E(g)$width <- 6 * module.module.rules$confidence + 0.5
    E(g)$color <- "gray30"
    E(g)$arrow.size <- 4
    
    # plot implication network
    
    layout(matrix(c(1, 2), 1, 2), c(2, 1), 1)
    par(mar=c(5, 4, 4, 1.8))
    
    image(matrix(spot.list$overview.mask, env$preferences$dim.1stLvlSom), col="gray90", axes=FALSE, main="Module implications (basket algorithm)", cex.main=1.6)
    title("all samples",line=0.2)
    box()
    
    par(new=TRUE, mar=c(3.9,2.7,2.7,0.5))
    plot(g, layout=t(sapply(spot.list$spots, function(x) { x$position })),
         vertex.label.color="black" ,vertex.label.cex=1, 
         xlim=c(-1.0,1.0), ylim=c(-1.0,1.0))
    
    par(mar=c(5, 1, 4, 2))
    plot(0, type="n", axes=FALSE, xlab="", ylab="", xlim=c(0,1), ylim=c(0,1), xaxs="i", yaxs="i")
    legend(0.08,0.95, c("Confidence","(group implication)","100%","50%","10%","","Confidence", "(module implication)","100%","50%","10%" ), 
           pch=c(NA,NA,1,1,1,rep(NA,6)), pt.cex=c(0,0,7.5,4,1,rep(0,6)), lty=c(rep(0,8),1,1,1), lwd=c(rep(0,8),7,4,1), bty="n", cex=1.2, y.intersp=2)
    box()
    
    
    for( gr in unique(env$group.labels) )
    {
      if( all( spotdata.binary[,which(env$group.labels==gr)] == 0 ) ) next
      
      module.module.rules <-
        do.call( rbind, lapply( rownames(spotdata.binary), function(mod.x.1)
        {
          do.call( rbind, lapply( setdiff( rownames(spotdata.binary), mod.x.1), function(mod.x.2)
          {
            count <- sum( colSums( spotdata.binary[c(mod.x.1,mod.x.2),which(env$group.labels==gr),drop=FALSE] ) == 2 )
            data.frame( lhs = mod.x.1,
                        rhs = mod.x.2,
                        count = count,
                        support = count / sum(env$group.labels==gr),
                        confidence = count / sum( spotdata.binary[mod.x.1,which(env$group.labels==gr)] == 1 ),
                        stringsAsFactors = FALSE)
          }) )
        } ) )
      module.module.rules <- module.module.rules[ which( module.module.rules$support>0.1 &
                                                           module.module.rules$confidence>0.1 ), ]
      module.module.rules <- module.module.rules[ order( module.module.rules$confidence, decreasing=TRUE), ]
      module.module.rules <- module.module.rules[ which( !duplicated( apply(cbind(module.module.rules$lhs,module.module.rules$rhs),1,function(x) paste(sort(x),collapse=" ") ) ) ),]
      
      
      if( nrow(module.module.rules) > 0 )   
      {
        g <- graph.empty( length(spot.list$spots), directed = FALSE )
        V(g)$name <- names(spot.list$spots)
        g <- add_edges( g, as.vector( rbind( match( module.module.rules$lhs, V(g)$name ), match( module.module.rules$rhs, V(g)$name ) ) ) )
        
        V(g)$label <- names(spot.list$spots)
        V(g)$size <- 10
        V(g)$color <- "gray"
        V(g)$color[ match( module.group.rules$lhs[which(module.group.rules$rhs==gr)], V(g)$name ) ] <- env$groupwise.group.colors[gr]   
        
        E(g)$width <- 6 * module.module.rules$confidence + 0.5
        E(g)$color <- env$groupwise.group.colors[gr]
        
        layout(matrix(c(1, 2), 1, 2), c(2, 1), 1)
        par(mar=c(5, 4, 4, 1.8))
        
        image(matrix(spot.list$overview.mask, env$preferences$dim.1stLvlSom), col="gray90", axes=FALSE, main="Module implications (basket algorithm)", cex.main=1.6)
        title(paste("Group samples:",gr),line=0.2,col.main=env$groupwise.group.colors[gr])
        box()
        
        par(new=TRUE, mar=c(3.9,2.7,2.7,0.5))
        plot(g, layout=t(sapply(spot.list$spots, function(x) { x$position })),
             vertex.label.color="black" ,vertex.label.cex=1, 
             xlim=c(-1.0,1.0), ylim=c(-1.0,1.0))
        
        par(mar=c(5, 1, 4, 2))
        plot(0, type="n", axes=FALSE, xlab="", ylab="", xlim=c(0,1), ylim=c(0,1), xaxs="i", yaxs="i")
        legend(0.08,0.95, c("Confidence","(module implication)","100%","50%","10%" ), 
               lty=c(rep(0,2),1,1,1), lwd=c(rep(0,2),7,4,1), bty="n", cex=1.2, y.intersp=2)
        box()
      }
    }
    
  }  
  

  #### Group Associations #### 
  
  if (length(unique(env$group.labels)) > 1)
  {
    layout( matrix(c(rep(1,4),2:13),nrow=4,byrow=TRUE), heights=c(0.15,1,1,1) )
    par(mar=c(0,0,0,0))
    
    frame()
    text(0.5,0.5,"Group associations of the modules",cex=2.6)
    
    i <- 1
    for( rowname in rownames(spotdata.binary) )
    { 
      if(i>1&&i%%13==0) { par(mar=c(0,0,0,0)); frame(); i <- i + 1 }
  
      count <- tapply(spotdata.binary[rowname,],env$group.labels,sum)[unique(env$group.labels)]
      percent <- count/table(env$group.labels)[unique(env$group.labels)]
      count <- count[which(count>0)]
      percent <- percent[which(percent>0)]
      
      if( length(percent) > 0 )
      { 
        par(mar=c(5,6,4,5))
        barplot( 100*percent, horiz=TRUE, main=paste(rowname," (",sum(count),")",sep=""), cex.main=1.6, xlim=c(0,100),
                 names.arg=paste(names(percent),"\n(",count,")",sep="" ), las=1, xlab="association (%)",
                 col=env$groupwise.group.colors[names(percent)] )
        i <- i + 1 
      } 
    }
  }  


  #### Group Implications ####
    
  if ( length(spot.list$spots) > 2 )
  {
    group.module.rules <-
      do.call( rbind, lapply( unique( env$group.labels ), function(gr.x)
      {
        samples.with.group <- names( which( env$group.labels == gr.x ) )
        do.call( rbind, lapply( rownames(spotdata.binary), function(mod.x)
        {
          count <- sum( spotdata.binary[mod.x,samples.with.group ] == 1 )
          data.frame( lhs = gr.x,
                      rhs = mod.x,
                      count = count,
                      support = count / ncol(env$indata),
                      confidence = count / length(samples.with.group),
                      stringsAsFactors = FALSE)
        }) )
      } ) )
    group.module.rules <- group.module.rules[ which( group.module.rules$support>0.01 &
                                                       group.module.rules$confidence>0.01 ), ]
    group.module.rules <- group.module.rules[ order( group.module.rules$confidence), ]
    
    
    layout( matrix(c(rep(1,4),2:13),nrow=4,byrow=TRUE), heights=c(0.15,1,1,1) )
    par(mar=c(0,0,0,0))
    
    frame()
    text(0.5,0.5,"Module implications of the groups (basket algorithm)",cex=2.6)
    
    i <- 1
    for( gr in unique(env$group.labels) )
    {
      if(i>1&&i%%13==0) { par(mar=c(0,0,0,0)); frame(); i <- i + 1 }
     
      r <- which( group.module.rules$lhs == gr )
      
      if( length(r) > 0 )
      {
        par(mar=c(5,6,4,5))
        barplot( 100*group.module.rules[r,]$confidence, horiz=TRUE, main=paste(gr," (",sum(env$group.labels==gr),")",sep=""), col.main=env$groupwise.group.colors[gr], cex.main=1.6, xlim=c(0,100),
                 names.arg=paste(group.module.rules[r,]$rhs,"\n(",group.module.rules[r,]$count,")",sep="" ), las=1, xlab="confidence" )
        title(main="implies",line=0,cex.main=0.8)
        i <- i + 1
      } 
    }
    
    
    module.group.rules <-
      do.call( rbind, lapply( rownames(spotdata.binary), function(mod.x)
      {
        samples.with.module <- names( which( spotdata.binary[mod.x,] == 1 ) )
        do.call( rbind, lapply( unique( env$group.labels[ samples.with.module ] ), function(gr.x)
        {
          count <- sum( env$group.labels[ samples.with.module ] == gr.x )
          data.frame( lhs = mod.x,
                      rhs = gr.x,
                      count = count,
                      support = count / ncol(env$indata),
                      confidence = count / length(samples.with.module),
                      stringsAsFactors = FALSE)
        }) )
      } ) )
    module.group.rules <- module.group.rules[ which( module.group.rules$support>0.01 &
                                                       module.group.rules$confidence>0.01 ), ]
    module.group.rules <- module.group.rules[ order( module.group.rules$confidence ), ]
    
    
    layout( matrix(c(rep(1,4),2:13),nrow=4,byrow=TRUE), heights=c(0.15,1,1,1) )  
    par(mar=c(0,0,0,0))
    
    frame()
    text(0.5,0.5,"Group implications of the modules (basket algorithm)",cex=2.6)
    
    i <- 1
    for( lhs in unique(rev(module.group.rules$lhs) ) )
    {
      if(i>1&&i%%13==0) { par(mar=c(0,0,0,0)); frame(); i <- i + 1 }
      
      r <- which( module.group.rules$lhs == lhs )
      n <- module.group.rules[r[1],]$count / module.group.rules[r[1],]$confidence
      
      if( length(r) > 0 )
      {
        par(mar=c(5,8,4,5))
        barplot( 100*module.group.rules[r,]$confidence, horiz=TRUE, main=paste(lhs," (",n,")",sep=""), cex.main=1.6, xlim=c(0,100),
                 names.arg=paste(module.group.rules[r,]$rhs,"\n(",module.group.rules[r,]$count,")",sep="" ), las=1, xlab="confidence",
                 col=env$groupwise.group.colors[module.group.rules[r,]$rhs] )
        title(main="implies",line=0,cex.main=0.8)
        i <- i + 1
      } 
    }
    
  }
    
  dev.off()
}
