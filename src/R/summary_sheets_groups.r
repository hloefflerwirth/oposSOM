pipeline.summarySheetsGroups <- function()
{
  group.metadata <- do.call(cbind, by(t(metadata), group.labels, colMeans))[,unique(group.labels)]
  group.metadata.sd <- do.call(cbind, by(t(metadata), group.labels, function(x) { apply(x,2,sd) }))[,unique(group.labels)]

  loglog.group.metadata <- apply(group.metadata, 2, function(x)
  {
    meta.sign <- sign(x)
    meta <- log10(abs(x))
    meta <- meta - min(meta, na.rm=TRUE)
    return(meta * meta.sign)
  })
  WAD.group.metadata <- apply(group.metadata ,2, function(x) { x * ((x - min(x)) / (max(x) - min(x))) })

  bleached.group.metadata <- group.metadata
  bleached.WAD.group.metadata <- WAD.group.metadata
  bleached.loglog.group.metadata <- loglog.group.metadata

  for (i in seq_along(unique(group.labels)))
  {
    pos.metagenes <- which(group.metadata[,i] >= 0)
    neg.metagenes <- which(group.metadata[,i] < 0)

    bleached.group.metadata[pos.metagenes,i] <-
      bleached.group.metadata[pos.metagenes,i] -
      pmin(bleached.group.metadata[pos.metagenes,i],
           apply(group.metadata[pos.metagenes,-i,drop=FALSE], 1, max))

    bleached.group.metadata[neg.metagenes,i] <-
      bleached.group.metadata[neg.metagenes,i] -
      pmax(bleached.group.metadata[neg.metagenes,i],
           apply(group.metadata[neg.metagenes,-i,drop=FALSE], 1, min))

    pos.metagenes <- which(WAD.group.metadata[,i] >= 0)
    neg.metagenes <- which(WAD.group.metadata[,i] < 0)

    bleached.WAD.group.metadata[pos.metagenes,i] <-
      bleached.WAD.group.metadata[pos.metagenes,i] -
      pmin(bleached.WAD.group.metadata[pos.metagenes,i],
           apply(WAD.group.metadata[pos.metagenes,-i,drop=FALSE], 1, max))

    bleached.WAD.group.metadata[neg.metagenes,i] <-
      bleached.WAD.group.metadata[neg.metagenes,i] -
      pmax(bleached.WAD.group.metadata[neg.metagenes,i],
           apply(WAD.group.metadata[neg.metagenes,-i,drop=FALSE], 1, min))

    pos.metagenes <- which(loglog.group.metadata[,i] >= 0)
    neg.metagenes <- which(loglog.group.metadata[,i] < 0)

    bleached.loglog.group.metadata[pos.metagenes,i] <-
      bleached.loglog.group.metadata[pos.metagenes,i] -
      pmin(bleached.loglog.group.metadata[pos.metagenes,i],
           apply(loglog.group.metadata[pos.metagenes,-i,drop=FALSE], 1, max))

    bleached.loglog.group.metadata[neg.metagenes,i] <-
      bleached.loglog.group.metadata[neg.metagenes,i] -
      pmax(bleached.loglog.group.metadata[neg.metagenes,i],
           apply(loglog.group.metadata[neg.metagenes,-i,drop=FALSE], 1, min))

  }

  filename <- file.path(paste(files.name, "- Results"),
                        "Summary Sheets - Groups",
                        "Expression Portraits Groups.pdf")

  util.info("Writing:", filename)
  pdf(filename, 29.7/2.54, 21/2.54)

  par(mfrow=c(5,6))
  par(mar=c(0.5,2.5,4.5,1.5))

  for (i in seq_along(unique(group.labels)))
  {
    plot(0, type="n", axes=FALSE, xlab="", ylab="", xlim=c(0,1))

    mtext(unique(group.labels)[i], side=3, line = 2, cex=1.5, at=-0.04, font=3,
          adj=0, col=groupwise.group.colors[i])

    par(new=TRUE)

    image(matrix(group.metadata[,i], preferences$dim.1stLvlSom, preferences$dim.1stLvlSom),
          axes=FALSE, col = colramp(1000))

    title(main="logFC", cex.main=1.5, line=0.5)
    box()

    image(matrix(bleached.group.metadata[,i],
                 preferences$dim.1stLvlSom,
                 preferences$dim.1stLvlSom),
          axes=FALSE, col = colramp(1000), zlim=range(bleached.group.metadata))

    title(main="group specific logFC", cex.main=1.3, line=0.5)
    box()

    image(matrix(WAD.group.metadata[,i],
                 preferences$dim.1stLvlSom,
                 preferences$dim.1stLvlSom),
          axes=FALSE, col = colramp(1000))

    title(main="WAD", cex.main=1.5, line=0.5)
    box()

    image(matrix(bleached.WAD.group.metadata[,i],
                 preferences$dim.1stLvlSom,
                 preferences$dim.1stLvlSom),
          axes=FALSE, col = colramp(1000), zlim=range(bleached.WAD.group.metadata))

    title(main="group specific WAD", cex.main=1.3, line=0.5)
    box()

    image(matrix(loglog.group.metadata[,i],
                 preferences$dim.1stLvlSom,
                 preferences$dim.1stLvlSom),
          axes=FALSE, col = colramp(1000))
    title(main="loglogFC", cex.main=1.5, line=0.5)
    box()

    image(matrix(bleached.loglog.group.metadata[,i],
                 preferences$dim.1stLvlSom,
                 preferences$dim.1stLvlSom),
          axes=FALSE, col = colramp(1000), zlim=range(bleached.loglog.group.metadata))

    title(main="group specific loglogFC", cex.main=1.3, line=0.5)
    box()
  }


  if (length(unique(group.labels)) <= 12)
  {
    # differential portraits
    l.matrix <- matrix(seq_len(length(unique(group.labels))^2), length(unique(group.labels)), byrow=TRUE)
    l.matrix <- cbind(rep(2, length(unique(group.labels))), l.matrix+2)
    l.matrix[1:(nrow(l.matrix)/2),1] <- 1
    layout(l.matrix, widths=c(5, rep(95/length(unique(group.labels)),length(unique(group.labels)))))

    par(mar=c(0,0,0,0))
    plot(0, type="n", axes=FALSE, xlab="", ylab="")
    text(1, 0, "differential expression: relative scale", srt=90, cex=2)
    plot(0, type="n", axes=FALSE, xlab="", ylab="")

    par(mar=c(0.5,2,2,0.5))
    zlim = c(0, 0)

    for (g1 in seq_along(unique(group.labels)))
    {
      for (g2 in seq_along(unique(group.labels)))
      {

        diff.metadata <- group.metadata[,g2] - group.metadata[,g1]

        if (g2 == 1)
        {
          col <- groupwise.group.colors[g1]
          attributes(col) <- NULL
          par("col.lab"=col, "mgp"=c(1,0,0))

          image(matrix(diff.metadata,
                       preferences$dim.1stLvlSom,
                       preferences$dim.1stLvlSom),
                axes=FALSE, col = colramp(1000), ylab=unique(group.labels)[g1])

          par("col.lab"="black", "mgp"=c(3,1,0))
        } else
        {
          image(matrix(diff.metadata,
                       preferences$dim.1stLvlSom,
                       preferences$dim.1stLvlSom),
                axes=FALSE, col = colramp(1000))
        }

        box()

        if (g1 == 1)
        {
          title(unique(group.labels)[g2], col.main=groupwise.group.colors[g2])
        }

        if (min(diff.metadata) < zlim[1])
        {
          zlim[1] <- min(diff.metadata)
        }

        if (max(diff.metadata) > zlim[2])
        {
          zlim[2] <- max(diff.metadata)
        }
      }
    }

    par(mar=c(0,0,0,0))
    plot(0, type="n", axes=FALSE, xlab="", ylab="")
    text(1, 0, "differential expression: absolute scale", srt=90, cex=2)
    plot(0, type="n", axes=FALSE, xlab="", ylab="")
    par(mar=c(0.5,2,2,0.5))

    for (g1 in seq_along(unique(group.labels)))
    {
      for (g2 in seq_along(unique(group.labels)))
      {

        diff.metadata <- group.metadata[,g2] - group.metadata[,g1]

        if (g2 == 1)
        {
          col <- groupwise.group.colors[g1]
          attributes(col) <- NULL
          par("col.lab"=col, "mgp"=c(1,0,0))

          image(matrix(diff.metadata,
                       preferences$dim.1stLvlSom,
                       preferences$dim.1stLvlSom),
                axes=FALSE, col=colramp(1000), zlim=zlim, ylab=unique(group.labels)[g1])

          par("col.lab"="black", "mgp"=c(3,1,0))
        } else
        {
          image(matrix(diff.metadata,
                       preferences$dim.1stLvlSom,
                       preferences$dim.1stLvlSom),
                axes=FALSE, col = colramp(1000), zlim=zlim)
        }

        box()

        if (g1 == 1)
        {
          title(unique(group.labels)[g2], col.main=groupwise.group.colors[g2])
        }
      }
    }

    par(mar=c(0,0,0,0))
    plot(0, type="n", axes=FALSE, xlab="", ylab="")
    text(1, 0, "significance: log10(p-value)", srt=90, cex=2)
    par(mar=c(8,2,8,1.8))

    image(x=1,y=seq(0,6,length.out=1000), z=rbind(1:1000),
          col=colorRampPalette(rep(c("blue3","green","yellow","orange","red","darkred"), each=3))(1000),
          axes=FALSE, xlab="", ylab="")

    axis(2, c(0:6), -c(0:6), las=2, cex.axis=0.75)
    box()

    par(mar=c(0.5,2,2,0.5))

    null.differences <- sample(group.metadata, 1000000, replace=TRUE) - sample(group.metadata, 1000000, replace=TRUE)
    null.culdensity <- ecdf(abs(null.differences))

    all.p.values <- c()

    for (g1 in seq_along(unique(group.labels)))
    {
      for (g2 in seq_along(unique(group.labels)))
      {
        diff.metadata <- group.metadata[,g2] - group.metadata[,g1]
        diff.pvalues <- 1 - null.culdensity(abs(diff.metadata))
        diff.pvalues[which(diff.pvalues<1e-6)] <- 1e-6
        all.p.values <- c(all.p.values, diff.pvalues)

        if (g2 == 1)
        {
          col <- groupwise.group.colors[g1]
          attributes(col) <- NULL
          par("col.lab"=col, "mgp"=c(1,0,0))

          image(matrix(log10(diff.pvalues), preferences$dim.1stLvlSom), axes=FALSE,
                col=rev(colorRampPalette(rep(c("blue3","green","yellow","orange","red","darkred"), each=1))(1000)),
                zlim=c(-6,0), ylab=unique(group.labels)[g1])

          par("col.lab"="black", "mgp"=c(3,1,0))
        } else
        {
          image(matrix(log10(diff.pvalues), preferences$dim.1stLvlSom), axes=FALSE,
                col=rev(colorRampPalette(rep(c("blue3","green","yellow","orange","red","darkred"), each=1))(1000)),
                zlim=c(-6,0))
        }

        box()

        if (g1 == 1)
        {
          title(unique(group.labels)[g2], col.main=groupwise.group.colors[g2])
        }
      }
    }

    par(mar=c(0,0,0,0))
    plot(0, type="n", axes=FALSE, xlab="", ylab="")
    text(1, 0, "significance: fdr", srt=90, cex=2)
    par(mar=c(8,2,8,1.8))

    image(x=1,y=seq(0,6,length.out=1000), z=rbind(1:1000),
          col=colorRampPalette(rep(c("blue3","green","yellow","orange","red","darkred"), each=3))(1000),
          axes=FALSE, xlab="", ylab="")

    axis(2, c(0:6), c(1,0.5,0.4,0.3,0.2,0.1,0), las=2, cex.axis=0.75)
    box()

    par(mar=c(0.5,2,2,0.5))

    fdrtool.result <- fdrtool(all.p.values, statistic="pvalue", plot=FALSE, verbose=FALSE)

    for (g1 in seq_along(unique(group.labels)))
    {
      for (g2 in seq_along(unique(group.labels)))
      {

        diff.metadata <- group.metadata[,g2] - group.metadata[,g1]
        diff.pvalues <- 1 - null.culdensity(abs(diff.metadata))
        diff.pvalues[which(diff.pvalues<1e-6)] <- 1e-6

        diff.fdr <- fdrtool.result$lfdr[match(diff.pvalues, fdrtool.result$pval)]
        diff.fdr[which(diff.fdr>0.6)] <- 0.6

        if (g2 == 1)
        {
          col <- groupwise.group.colors[g1]
          attributes(col) <- NULL
          par("col.lab"=col, "mgp"=c(1,0,0))

          image(matrix(diff.fdr, preferences$dim.1stLvlSom), axes=FALSE,
                col=rev(colorRampPalette(rep(c("blue3","green","yellow","orange","red","darkred"), each=1))(1000)),
                zlim=c(0,0.6), ylab=unique(group.labels)[g1])

          par("col.lab"="black", "mgp"=c(3,1,0))
        } else
        {
          image(matrix(diff.fdr, preferences$dim.1stLvlSom), axes=FALSE,
                col=rev(colorRampPalette(rep(c("blue3","green","yellow","orange","red","darkred"), each=1))(1000)),
                zlim=c(0,0.6))
        }

        box()

        if (g1==1)
        {
          title(unique(group.labels)[g2], col.main=groupwise.group.colors[g2])
        }
      }
    }
  }

  dev.off()


  ###### Group stability scores
  filename <- file.path(paste(files.name, "- Results"),
                        "Summary Sheets - Groups",
                        "Group Assignment.pdf")

  util.info("Writing:", filename)
  pdf(filename, 21/2.54, 29.7/2.54)


  S <- tapply( group.silhouette.coef, group.labels, sort, decreasing=TRUE, simplify=FALSE )[unique(group.labels)]
  names(S) <- NA
  S <- unlist(S)
  names(S) <- sub( paste("^NA.",sep=""), "", names(S) )
  
  
  PCM <- cor( metadata )
  diag(PCM) <- NA
   
  group.correlations <- sapply( seq(ncol(metadata)), function(i)
  {
    mean.group.correlations <- tapply( PCM[,i], group.labels, mean, na.rm=TRUE )[unique(group.labels)]
    
    return(  mean.group.correlations )
  } )
  colnames(group.correlations) <- colnames(indata)
  group.correlations[which(is.nan(group.correlations))] <- 0
    
  
  layout(matrix(c(0,1,2,0),ncol=1))
  #par(mfrow=c(2,1))
  par(mar=c(5,3,3,2))
  
  b<-barplot( S, col=group.colors[names(S)], main="Correlation Silhouette", names.arg=if(ncol(indata)<80) names(S) else rep("",length(S)), las=2, cex.main=1, cex.lab=1, cex.axis=0.8, cex.names=0.6, border = ifelse(ncol(indata)<80,"black",NA), xpd=FALSE, ylim=c(-.25,1) )  
  mtext("S",2,line=1.9,cex=0.8)
  abline( h=c(0,0.25,0.5,0.75), lty=2, col="gray80" )
  title( main= bquote("<" ~ s ~ "> = " ~ .(round(mean(S),2))), line=0.5, cex.main=1 )
  box()
  points( b, rep(-0.2,ncol(indata)), pch=15, cex=1, col=groupwise.group.colors[apply( group.correlations[,names(S)], 2, which.max )] )
  
  mean.boxes <- by( S, group.labels, c )[ unique( group.labels ) ]
  mean.mean.S <- sapply( mean.boxes, mean )
  
  par(mar=c(5,3,0,2))
  boxplot( mean.boxes, col=groupwise.group.colors, las=2, main="", cex.main=1, cex.axis=0.8, xaxt="n", ylim=c(-.25,1) )
  mtext("S",2,line=1.9,cex=0.8)
  abline( h=c(0,0.25,0.5,0.75), lty=2, col="gray80" )
  axis( 1, 1:length(groupwise.group.colors), paste( unique(group.labels), "\n<s> =", round(mean.mean.S,2) ), las=2, cex.axis=0.8 )

 
  
  par(mfrow=c(1,1), mar=c(5,15,1,1) )
  
  for( gr in unique(group.labels) )
  {
    samples <- names(which(group.labels==gr))
    samples.o <- rev( names(S[which(names(S)%in%samples)]) )
    
    image( group.correlations[,samples.o,drop=F], col=colorRampPalette(c("gray90","orange","red4"))(1000) , zlim=c(0,1), axes=FALSE )
    box()	
    
    dummy<-sapply(1:length(unique(group.labels)), function(i)
    {	
      axis(1, seq(0,1,length.out=length(unique(group.labels)))[i], unique(group.labels)[i], col.axis=groupwise.group.colors[i], las=2, cex.axis=0.6)
    }	)
    
    for( i in seq(samples.o) )
    {
      axis(2, seq(0,1,length.out=length(samples.o))[i], samples.o[i], las=2, col.axis=groupwise.group.colors[gr], line=10, tick=FALSE, cex.axis=0.6 )      
      axis(2, seq(0,1,length.out=length(samples.o))[i], bquote("<" ~ r[.(gr)] ~ "> = " ~ .(round(group.correlations[gr,samples.o[i]],2)) ) , las=2, col.axis=groupwise.group.colors[gr], line=5, tick=FALSE, cex.axis=0.6 )
            
      second.corr.group <- sort( group.correlations[,samples.o[i]], decreasing=TRUE )      
      second.corr.group <- names(second.corr.group)[ which(names(second.corr.group)!=gr)[1] ]
      axis(2, seq(0,1,length.out=length(samples))[i], bquote("<" ~ r[.(second.corr.group)] ~ "> = " ~ .(round(group.correlations[second.corr.group,samples.o[i]],2)) ) , las=2, col.axis=groupwise.group.colors[second.corr.group], line=0, tick=FALSE, cex.axis=0.6 )					
    }
    
  }
  
  dev.off()

  
  
  

  ### Group clustering
  .count <- 0
  .topAbsis <- 0
  .heights <- 0

  absi <- function(hc, level=length(hc$height), init=TRUE)
  {
    if (init)
    {
      .count <<- 0
      .topAbsis <<- NULL
      .heights <<- NULL
    }

    if (level<0)
    {
      .count <<- .count + 1
      return(.count)
    }

    node <- hc$merge[level,]
    le <- absi(hc, node[1], init=FALSE)
    ri <- absi(hc, node[2], init=FALSE)
    mid <- (le+ri)/2
    .topAbsis <<- c(.topAbsis, mid)
    .heights <<- c(.heights, hc$height[level])
    invisible(mid)
  }

  filename <- file.path(paste(files.name, "- Results"),
                        "Summary Sheets - Groups",
                        "Group Clustering.pdf")

  util.info("Writing:", filename)
  pdf(filename , 29.7/2.54, 21/2.54)

  for (i in seq_along(unique(group.labels)))
  {
    group.member <- which(group.labels==unique(group.labels)[i])

    if (length(group.member) >= 5)
    {
      hc <- hclust(dist(t(metadata[,group.member])))#, method="average")
      hc$labels <- names(group.member)

      absi(hc)
      o <- order(.heights)
      coords.h <- .heights[o]
      coords.x <- .topAbsis[o]

      merges <- list()
      old.cluster.size <- rep(1,length(group.member))

      for (hi in c(seq_along(hc$height)))
      {
        new.clusters <- cutree(hc, h=hc$height[hi])
        new.cluster.size <- table(new.clusters)[new.clusters]
        merges[[hi]] <- which(old.cluster.size != new.cluster.size)
        old.cluster.size <- new.cluster.size
      }

      equal.merges <- which(sapply(merges,length) == 0)

      for (em.i in equal.merges)
      {
        root.branch <- max(which(setdiff(c(seq_along(merges)), equal.merges) < em.i))
        merges[[em.i]] <- merges[[root.branch]]
      }

      top.level <- c(0)

      for (hi in seq(2, length(hc$height)))
      {
        mem <- rev(merges)[[hi]]
        last.level <- 0

        for (hi2 in 1:(hi-1))
        {
          if (any(mem %in% rev(merges)[[hi2]]))
          {
            last.level <- hi2
          }
        }

        top.level[hi] <- last.level
      }

      top.level <- rev(top.level)
      diff.metadata <- list(rowMeans(metadata[, group.member]))

      for (hi in seq(2, length(hc$height)))
      {
        mem <- rev(merges)[[hi]]
        tl <- rev(top.level)[hi]

        diff.metadata[[hi]] <- rowMeans(metadata[, group.member[mem]])
        diff.metadata[[hi]][which(is.na(diff.metadata[[tl]]))] <- NA
        diff.metadata[[hi]][which(diff.metadata[[tl]] > quantile(diff.metadata[[tl]],0.9, na.rm=TRUE))] <- NA
        diff.metadata[[hi]][which(diff.metadata[[tl]] < quantile(diff.metadata[[tl]],0.1, na.rm=TRUE))] <- NA
      }

      diff.metadata <- rev(diff.metadata)

      if (length(group.member) >= 20)
      {
        cex.branch.portraits <- c(0.6, 2 * 0.6 * max(hc$height) / length(group.member))
        cex.sample.portraits <- c(0.4, 2 * 0.4 * max(hc$height) / length(group.member))
        y.sample.portraits <- -1
      } else
      {
        c1 <- 0.1 + (0.4 * (length(group.member)-5)) / 15
        c2 <- 0.1 + (0.3 * (length(group.member)-5)) / 15
        cex.branch.portraits <- c(c1, 2 * c1 * max(hc$height) / length(group.member))
        cex.sample.portraits <- c(c2, 2 * c2 * max(hc$height) / length(group.member))
        y.sample.portraits <- 0
      }

      plot(hc, main=unique(group.labels)[i], col.main=unique(group.colors)[i], xlab="", sub="")
      mtext("mean branch portraits")

      for (ii in (if (length(group.member)<80) 1 else ceiling(length(merges)/3)):length(merges))
      {
        m <- matrix(rowMeans(metadata[, group.member[merges[[ii]]]]), preferences$dim.1stLvlSom, preferences$dim.1stLvlSom)

        if (max(m) - min(m) != 0)
        {
          m <- 1 + (m - min(m)) / (max(m) - min(m)) * 999
        }

        m <- cbind(apply(m, 1, function(x){x}))[nrow(m):1,]
        x <- pixmapIndexed(m , col = colramp(1000), cellres=10)

        addlogo(x,
                coords.x[ii]+cex.branch.portraits[1]*c(-1,1),
                coords.h[ii]+cex.branch.portraits[2]*c(-1,1))

        rect(coords.x[ii] - cex.branch.portraits[1],
             coords.h[ii] - cex.branch.portraits[2],
             coords.x[ii] + cex.branch.portraits[1],
             coords.h[ii] + cex.branch.portraits[2])
      }

      if (length(group.member) < 80)
      {
        for (ii in seq_along(group.member))
        {
          m <- matrix(metadata[, group.member[hc$order[ii]]],
                      preferences$dim.1stLvlSom, preferences$dim.1stLvlSom)

          if (max(m) - min(m) != 0)
          {
            m <- 1 + (m - min(m)) / (max(m) - min(m)) * 999
          }

          m <- cbind(apply(m, 1, function(x){x}))[nrow(m):1,]
          x <- pixmapIndexed(m , col = colramp(1000), cellres=10)

          addlogo(x,
                  ii+cex.sample.portraits[1]*c(-1,1),
                  y.sample.portraits+cex.sample.portraits[2]*c(-1,1))
        }
      }

      plot(hc, main=unique(group.labels)[i], col.main=unique(group.colors)[i], xlab="", sub="")
      mtext("mean branch portraits; difference to mean group portrait")

      for (ii in (if (length(group.member)<80) 1 else ceiling(length(merges)/3)):length(merges))
      {
        m <- matrix(rowMeans(metadata[, group.member[merges[[ii]]]]) - rowMeans(metadata[, group.member]),
                    preferences$dim.1stLvlSom, preferences$dim.1stLvlSom)

        if (max(m,na.rm=TRUE) - min(m,na.rm=TRUE) != 0)
        {
          m <- 1 + (m - min(m,na.rm=TRUE)) / (max(m,na.rm=TRUE) - min(m,na.rm=TRUE)) * 999
        }

        m <- cbind(apply(m, 1, function(x){x}))[nrow(m):1,]
        m[which(is.na(m))] <- 0
        x <- pixmapIndexed(m , col = c("gray90", colramp(1000)), cellres=10)

        addlogo(x,
                coords.x[ii]+cex.branch.portraits[1]*c(-1,1),
                coords.h[ii]+cex.branch.portraits[2]*c(-1,1))

        rect(coords.x[ii] - cex.branch.portraits[1],
             coords.h[ii] - cex.branch.portraits[2],
             coords.x[ii] + cex.branch.portraits[1],
             coords.h[ii] + cex.branch.portraits[2])
      }

      if (length(group.member) < 80)
      {
        for (ii in seq_along(group.member))
        {
          m <- matrix(metadata[, group.member[hc$order[ii]]],
                      preferences$dim.1stLvlSom, preferences$dim.1stLvlSom)

          if (max(m) - min(m) != 0)
          {
            m <- 1 + (m - min(m)) / (max(m) - min(m)) * 999
          }

          m <- cbind(apply(m, 1, function(x){x}))[nrow(m):1,]
          x <- pixmapIndexed(m , col = colramp(1000), cellres=10)

          addlogo(x,
                  ii+cex.sample.portraits[1]*c(-1,1),
                  y.sample.portraits+cex.sample.portraits[2]*c(-1,1))
        }
      }

    }
  }

  dev.off()
}
