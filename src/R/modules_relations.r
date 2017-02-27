modules.relations <- function(spot.list, main, path)
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
    
    image(matrix(spot.list$overview.mask, preferences$dim.1stLvlSom), col="gray90", axes=FALSE, main="Module correlations (WTO algorithm)", cex.main=1.6)
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
  
  if ( length(spot.list$spots) > 2 )
  {
    trans.groups <- as(as.data.frame(group.labels), "transactions")
    trans.modules <- as(t(spotdata.binary), "transactions")
    trans.merge <- merge(trans.groups, trans.modules)
    
    suppressWarnings({
      rules <- apriori(trans.merge, parameter = list(supp = 0.01, minlen = 2, maxlen = 2,
                                                   conf = 0.1, target = "rules"),
                     control = list(verbose=FALSE))  })
    rules <- as(rules,"data.frame")
    rules <- rules[ order(rules$confidence,decreasing = TRUE), ]  
    rules$lhs <- sapply( strsplit(as.character(rules$rules)," => ", fixed=TRUE), function(x) gsub("[{}]","",x) )[1,]
    rules$rhs <- sapply( strsplit(as.character(rules$rules)," => ", fixed=TRUE), function(x) gsub("[{}]","",x) )[2,]
    
    rules.modules <- rules[ -unique( c( grep("group.labels=",rules$lhs), grep("group.labels=",rules$rhs) )),]
    rules.modules <- rules.modules[ which( !duplicated( apply(cbind(rules.modules$lhs,rules.modules$rhs),1,function(x) paste(sort(x),collapse=" ") ) ) ),]
    rules.groups <- rules[ grep("group.labels=",rules$rhs) ,]
    rules.groups <- rules.groups[ !duplicated( rules.groups$lhs ) ,]
    rules.groups$rhs <- sub( "group.labels=","",rules.groups$rhs )
    
    g <- graph.empty( length(spot.list$spots), directed = FALSE )
    V(g)$name <- names(spot.list$spots)
    g <- add_edges( g, as.vector( rbind( match( rules.modules$lhs, V(g)$name ), match( rules.modules$rhs, V(g)$name ) ) ) )
    V(g)$label <- names(spot.list$spots)
    V(g)$size <- 5
    V(g)$size[ match( rules.groups$lhs, V(g)$name ) ] <- 15 * rules.groups$confidence + 3
    V(g)$color <- "gray"
    V(g)$color[ match( rules.groups$lhs, V(g)$name ) ] <- groupwise.group.colors[rules.groups$rhs]
    E(g)$width <- 6 * rules.modules$confidence + 0.5
    E(g)$color <- "gray30"
    E(g)$arrow.size <- 4
    
    # plot implication network
    
    layout(matrix(c(1, 2), 1, 2), c(2, 1), 1)
    par(mar=c(5, 4, 4, 1.8))
    
    image(matrix(spot.list$overview.mask, preferences$dim.1stLvlSom), col="gray90", axes=FALSE, main="Module implications (basket algorithm)", cex.main=1.6)
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
    
    for( gr in unique(group.labels) )
    {
      if( all( spotdata.binary[,which(group.labels==gr)] == 0 ) ) next
            
      trans.modules <- as(t(spotdata.binary[,which(group.labels==gr)]), "transactions")
      rules <- suppressWarnings({ apriori(trans.modules, parameter = list(supp = 0.1, minlen = 2, maxlen = 2,
                                                                          conf = 0.1, target = "rules"),
                                          control = list(verbose=FALSE))  })
      if( length(rules) > 0 )
      {
        rules <- as(rules,"data.frame")
        rules <- rules[ order(rules$confidence,decreasing = TRUE), ]  
        rules$lhs <- sapply( strsplit(as.character(rules$rules)," => ", fixed=TRUE), function(x) gsub("[{}]","",x) )[1,]
        rules$rhs <- sapply( strsplit(as.character(rules$rules)," => ", fixed=TRUE), function(x) gsub("[{}]","",x) )[2,]
        
        rules.modules <- rules[ which( !duplicated( apply(cbind(rules$lhs,rules$rhs),1,function(x) paste(sort(x),collapse=" ") ) ) ),]
        
        
        g <- graph.empty( length(spot.list$spots), directed = FALSE )
        V(g)$name <- names(spot.list$spots)
        g <- add_edges( g, as.vector( rbind( match( rules.modules$lhs, V(g)$name ), match( rules.modules$rhs, V(g)$name ) ) ) )
        
        V(g)$label <- names(spot.list$spots)
        V(g)$size <- 10
        V(g)$color <- "gray"
        V(g)$color[ match( rules.groups$lhs[which(rules.groups$rhs==gr)], V(g)$name ) ] <- groupwise.group.colors[gr]   
        
        E(g)$width <- 6 * rules.modules$confidence + 0.5
        E(g)$color <- groupwise.group.colors[gr]
        
        layout(matrix(c(1, 2), 1, 2), c(2, 1), 1)
        par(mar=c(5, 4, 4, 1.8))
        
        image(matrix(spot.list$overview.mask, preferences$dim.1stLvlSom), col="gray90", axes=FALSE, main="Module implications (basket algorithm)", cex.main=1.6)
        title(paste("Group samples:",gr),line=0.2,col.main=groupwise.group.colors[gr])
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
  
  if (length(unique(group.labels)) > 1)
  {
    layout( matrix(c(rep(1,4),2:13),nrow=4,byrow=TRUE), heights=c(0.15,1,1,1) )
    par(mar=c(0,0,0,0))
    
    frame()
    text(0.5,0.5,"Group associations of the modules",cex=2.6)
    
    i <- 1
    for( rowname in rownames(spotdata.binary) )
    { 
      if(i>1&&i%%13==0) { par(mar=c(0,0,0,0)); frame(); i <- i + 1 }
  
      count <- tapply(spotdata.binary[rowname,],group.labels,sum)[unique(group.labels)]
      percent <- count/table(group.labels)[unique(group.labels)]
      count <- count[which(count>0)]
      percent <- percent[which(percent>0)]
      
      if( length(percent) > 0 )
      { 
        par(mar=c(5,6,4,5))
        barplot( 100*percent, horiz=TRUE, main=paste(rowname," (",sum(count),")",sep=""), cex.main=1.6, xlim=c(0,100),
                 names.arg=paste(names(percent),"\n(",count,")",sep="" ), las=1, xlab="association (%)",
                 col=groupwise.group.colors[names(percent)] )
        i <- i + 1 
      } 
    }
  }  


  #### Group Implications ####
    
  if ( length(spot.list$spots) > 2 )
  {
    trans.groups <- as(as.data.frame(group.labels), "transactions")
    trans.modules <- as(t(spotdata.binary), "transactions")
    trans.merge <- merge(trans.groups, trans.modules)
    
    suppressWarnings({
      rules <- apriori(trans.merge, parameter = list(supp = 0.01, minlen = 2, maxlen = nrow(spotdata.binary),
                                                   conf = 0.01, target = "rules"),
                     control = list(verbose=FALSE))  })
    rules <- as(rules,"data.frame")
    rules$n <- rules$support * length(trans.merge)
    
    rules$lhs <- sapply( strsplit(as.character(rules$rules)," => ", fixed=TRUE), function(x) gsub("[{}]","",x) )[1,]
    rules$rhs <- sapply( strsplit(as.character(rules$rules)," => ", fixed=TRUE), function(x) gsub("[{}]","",x) )[2,]
    
    rules <- rules[ order(rules$confidence), ]  
    
    
    layout( matrix(c(rep(1,4),2:13),nrow=4,byrow=TRUE), heights=c(0.15,1,1,1) )
    par(mar=c(0,0,0,0))
    
    frame()
    text(0.5,0.5,"Module implications of the groups (basket algorithm)",cex=2.6)
    
    i <- 1
    for( gr in unique(group.labels) )
    {
      if(i>1&&i%%13==0) { par(mar=c(0,0,0,0)); frame(); i <- i + 1 }
      
      r <- which( rules$lhs == paste("group.labels=",gr,sep="") )
      
      if( length(r) > 0 )
      {
        par(mar=c(5,6,4,5))
        barplot( 100*rules[r,]$confidence, horiz=TRUE, main=paste(gr," (",sum(group.labels==gr),")",sep=""), col.main=groupwise.group.colors[gr], cex.main=1.6, xlim=c(0,100),
                 names.arg=paste(rules[r,]$rhs,"\n(",rules[r,]$n,")",sep="" ), las=1, xlab="confidence" )
        title(main="implies",line=0,cex.main=0.8)
        i <- i + 1
      } 
    }
    
    rules.help <- rules[ which(rules$rhs %in% paste("group.labels=",unique(group.labels),sep="") ),]
    rules.help$rhs <- sub( "group.labels=","",rules.help$rhs )
    
    layout( matrix(c(rep(1,4),2:13),nrow=4,byrow=TRUE), heights=c(0.15,1,1,1) )  
    par(mar=c(0,0,0,0))
    
    frame()
    text(0.5,0.5,"Group implications of the modules (basket algorithm)",cex=2.6)
    
    i <- 1
    for( lhs in unique(rev(rules.help$lhs) ) )
    {
      if(i>1&&i%%13==0) { par(mar=c(0,0,0,0)); frame(); i <- i + 1 }
      
      r <- which( rules.help$lhs == lhs )
      n <- rules.help[r[1],]$n/ rules.help[r[1],]$confidence
      
      if( length(r) > 0 )
      {
        par(mar=c(5,8,4,5))
        barplot( 100*rules.help[r,]$confidence, horiz=TRUE, main=paste(lhs," (",n,")",sep=""), cex.main=1.6, xlim=c(0,100),
                 names.arg=paste(rules.help[r,]$rhs,"\n(",rules.help[r,]$n,")",sep="" ), las=1, xlab="confidence",
                 col=groupwise.group.colors[rules.help[r,]$rhs] )
        title(main="implies",line=0,cex.main=0.8)
        i <- i + 1
      } 
    }
    
  }
    
  dev.off()
}
