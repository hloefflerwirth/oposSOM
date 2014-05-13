pipeline.supportingMaps <- function()
{
  dirname <- paste(files.name, "- Results/Supporting Maps&Profiles")
  dir.create(dirname, showWarnings=F)

  filename <- file.path(dirname, "Supporting Maps.pdf")
  util.info("Writing:", filename)

  pdf(filename, 21/2.54, 21/2.54)

  ### Population Map ###
  par(mfrow = c(1, 1))
  par(mar=c(5, 6, 4, 5))
  m <- log10(som.result$code.sum[,"nobs"])
  m[which(is.infinite(m))] <- NA

  image(matrix(m, preferences$dim.1stLvlSom, preferences$dim.1stLvlSom),
        axes=F, col=colramp(1000), main="Population Map", cex.main=2.5)

  mtext("log ( # genes in metagene )", side=1, line=1, cex=1.4)
  box()

  par(new=T, mar=c(1, 0, 0, 0))
  layout(matrix(c(0, 0, 0, 0, 1, 0, 0, 0, 0), 3, 3), c(1, 0.05, 0.02), c(0.14, 0.3, 1))
  image(matrix(1:100, 1, 100), col = colramp(1000), axes=F)

  axis(2, c(1, max(som.result$code.sum[,"nobs"])),
       at=c(0, 1), las=2, tick=F, pos=-0.5, cex.axis=1.4)

  box()

  ### Variance Map ###
  par(mfrow = c(1, 1))
  par(mar=c(5, 6, 4, 5))

  image(matrix(log10(apply(metadata, 1, var)), preferences$dim.1stLvlSom, preferences$dim.1stLvlSom),
        axes=F, col=colramp(1000), main="Metagene Variance Map", cex.main=2.5)

  mtext("log ( metagene variance )", side=1, line=1, cex=1.4)
  box()

  par(new=T, mar=c(1, 0, 0, 0))
  layout(matrix(c(0, 0, 0, 0, 1, 0, 0, 0, 0), 3, 3), c(1, 0.05, 0.02), c(0.14, 0.3, 1))
  image(matrix(1:100, 1, 100), col=colramp(1000), axes=F)

  axis(2, c(round(min(apply(metadata, 1, var), na.rm=T),2),
            round(max(apply(metadata, 1, var),na.rm=T),2)),
       at=c(0, 1), las=2, tick=F, pos=-0.5, cex.axis=1.4)

  box()

  ### Significance Map ###
  par(mfrow = c(1, 1))
  par(mar=c(5, 6, 4, 5))
  suppressWarnings({ p <- -apply(p.m, 1, mean) })

  image(matrix(p, preferences$dim.1stLvlSom, preferences$dim.1stLvlSom),
        axes=F, main="Metagene Significance Map", cex.main=2.5, zlim=c(-1,0),
        col=colorRampPalette(c("blue4", "blue4", "blue3", "blue3", "blue2",
                               "blue2", "blue1", "lightblue", "darkgreen",
                               "#008B00", "green3", "green", "yellow", "gold",
                               "orange", "red", "darkred"))(1000))

  mtext(expression(paste("<",p[k],">")), side=1, line=1, cex=1.4)
  box()

  par(new=T, mar=c(1, 0, 0, 0))
  layout(matrix(c(0, 0, 0, 0, 1, 0, 0, 0, 0), 3, 3), c(1, 0.05, 0.02), c(0.14, 1.15, 0.15))

  image(matrix(1:10, 1, 10), axes=F,
        col=rev(colorRampPalette(c("darkgreen", "#008B00", "green3", "green",
                                   "yellow", "gold", "orange", "red", "darkred"))(1000)))

  axis(2, c("0.0", "0.05", "0.1", "0.2", "0.3", "0.4", "0.5"),
       at=c(-0.05, 0.06, 0.17, 0.39, 0.61, 0.83, 1.05), las=2, tick=F, pos=-0.5, cex.axis=1.4)

  box()

  par(mfrow = c(1, 1))
  par(mar=c(5, 6, 4, 5))
  suppressWarnings({ p <- apply(p.m, 1, min) })

  image(matrix(-log(p), preferences$dim.1stLvlSom, preferences$dim.1stLvlSom),
        axes=F, main="Metagene Significance Map", cex.main=2.5,
        col=colorRampPalette(c("blue4", "blue4", "blue3", "blue3", "blue2",
                               "blue2", "blue1", "lightblue", "darkgreen",
                               "#008B00", "green3", "green", "yellow", "gold",
                               "orange", "red", "darkred"))(1000))

  mtext(expression(paste(min[m],"( log ",p["k,m"]," )")), side=1, line=1, cex=1.4)
  box()

  par(new=T, mar=c(1, 0, 0, 0))
  layout(matrix(c(0, 0, 0, 0, 1, 0, 0, 0, 0), 3, 3), c(1, 0.05, 0.02), c(0.14, 0.3, 1))
  image(matrix(1:100, 1, 100), col=rev(colramp(1000)), axes=F)
  axis(2, at=c(0,1), c(round(min(log(p),na.rm=T)), 0), las=2, tick=F, pos=-0.5, cex.axis=1.4)
  box()


  ### Entropy Map ###
  q25 <- quantile(metadata, 0.25)
  q75 <- quantile(metadata, 0.75)

  p.metadata <- apply(metadata, 1, function(x)
  {
    hist(x, breaks=c(min(x), q25, q75, max(x)), plot=F)$counts / ncol(indata)
  })

  H <- apply(p.metadata, 2, function(p)
  {
    p <- p[p != 0]
    -sum(p * log2(p))
  })

  par(mfrow = c(1, 1))
  par(mar=c(5, 6, 4, 5))

  image(matrix(H, preferences$dim.1stLvlSom, preferences$dim.1stLvlSom),
        axes=F, col=colramp(1000), main="Standard Metagene Entropy Map", cex.main=2.5)

  mtext(expression(h[k]), side=1, line=1, cex=1.4)
  box()

  par(new=T, mar=c(1, 0, 0, 0))
  layout(matrix(c(0, 0, 0, 0, 1, 0, 0, 0, 0), 3, 3), c(1, 0.05, 0.02), c(0.14, 0.3, 1))
  image(matrix(1:100, 1, 100), col = colramp(1000), axes=F)

  axis(2, c(round(min(H,na.rm=T),2),round(max(H,na.rm=T),2)),
       at=c(0, 1), las=2, tick=F, pos=-0.5, cex.axis=1.4)

  box()

  ### Covariance Map ###
  errors <- c()

  for (i in 1:preferences$dim.1stLvlSom^2)
  {
    genes <- names(which(som.nodes == i))
    mean.cor <- 0

    for (ii in genes)
    {
      suppressWarnings({
        mean.cor <- mean.cor + cor(metadata[i,], as.numeric(indata[ii,]))
      })
    }

    errors[i] <- mean.cor / length(genes)
  }

  par(mfrow = c(1, 1))
  par(mar=c(5, 6, 4, 5))

  image(matrix(errors, preferences$dim.1stLvlSom, preferences$dim.1stLvlSom), axes=F,
        col=colramp(1000), main="Gene-Metagene Covariance Map", cex.main=2.5)

  mtext("correlation genes - metagene", side=1, line=1, cex=1.4)
  box()

  par(new=T, mar=c(1, 0, 0, 0))
  layout(matrix(c(0, 0, 0, 0, 1, 0, 0, 0, 0), 3, 3), c(1, 0.05, 0.02), c(0.14, 0.3, 1))
  image(matrix(1:100, 1, 100), col = colramp(1000), axes=F)

  axis(2, c(round(min(errors,na.rm=T),2),round(max(errors,na.rm=T),2)),
       at=c(0, 1), las=2, tick=F, pos=-0.5, cex.axis=1.4)

  box()

  ### Deviation Map ###
  errors <- rep(NA, preferences$dim.1stLvlSom^2)

  for (i in 1:preferences$dim.1stLvlSom^2)
  {
    genes <- names(which(som.nodes == i))

    e <- apply(indata[genes, ,drop=F], 1, function(x)
    {
      1 / (ncol(indata)-1) * sum((x - metadata[i,])^2)
    })

    errors[i] <- sqrt(mean(e))
  }

  par(mfrow = c(1, 1))
  par(mar=c(5, 6, 4, 5))

  image(matrix(errors, preferences$dim.1stLvlSom, preferences$dim.1stLvlSom),
        axes=F, col=colramp(1000), main="Deviation Map", cex.main=2.5)

  mtext("deviation genes - metagene", side=1, line=1, cex=1.4)
  box()

  par(new=T, mar=c(1, 0, 0, 0))
  layout(matrix(c(0, 0, 0, 0, 1, 0, 0, 0, 0), 3, 3), c(1, 0.05, 0.02), c(0.14, 0.3, 1))
  image(matrix(1:100, 1, 100), col=colramp(1000), axes=F)

  axis(2, c(round(min(errors,na.rm=T),2), round(max(errors,na.rm=T),2)),
       at=c(0,1), las=2, tick=F, pos=-0.5, cex.axis=1.4)

  box()

  ### Distance Map ###
  uh <- rep(NA, preferences$dim.1stLvlSom^2)

  for (i in 1:preferences$dim.1stLvlSom^2)
  {
    pos.x <- som.result$code.sum[i,1] + 1
    pos.y <- som.result$code.sum[i,2] + 1

    uh[i] <- mean(sapply(get.neighbors(pos.x, pos.y, preferences$dim.1stLvlSom), function(x)
    {
      x <- x - 1
      ij <- (x[1] + 1) + x[2] * preferences$dim.1stLvlSom
      sqrt(sum((metadata[ij,] - metadata[i,])^2))
    }))
  }

  par(mfrow = c(1, 1))
  par(mar=c(5, 6, 4, 5))

  image(matrix(-log10((uh)), preferences$dim.1stLvlSom, preferences$dim.1stLvlSom),
        axes=F, col=colorRampPalette(c("white","black"))(1000),
        main="Distance Map", cex.main=2.5)

  mtext("deviation of adjacent metagenes", side=1, line=1, cex=1.4)
  box()

  par(new=T, mar=c(1, 0, 0, 0))
  layout(matrix(c(0, 0, 0, 0, 1, 0, 0, 0, 0), 3, 3), c(1, 0.05, 0.02), c(0.14, 0.3, 1))
  image(matrix(1:100, 1, 100), col = colorRampPalette(c("white", "black"))(1000), axes=F)

  axis(2, c(round(min(uh,na.rm=T),2), round(max(uh,na.rm=T),2)),
       at=c(0,1), las=2, tick=F, pos=-0.5, cex.axis=1.4)

  box()

  ### Plot SOM 10x10 grid ###
  co.dim <- 10

  if (preferences$dim.1stLvlSom != co.dim)
  {
    co.so <- som.init(indata, xdim=10, ydim=10)

    # Rotate First lvl SOMs
    if (preferences$rotate.SOM.portraits > 0)
    {
      for (i in 1:preferences$rotate.SOM.portraits)
      {
        o <- matrix(c(1:(100)), 10, 10, byrow=T)
        o <- o[rev(1:10) ,]
        co.so <- co.so[as.vector(o),]
      }
    }

    # Train SOM
    co.so <- som.train(indata, co.so, xdim=10, ydim=10, alpha=0.05, radius=10,
                       rlen=nrow(indata)*2, inv.alp.c=nrow(indata)*2/100)

    co.so <- som.train(indata, co.so$code, xdim=10, ydim=10, alpha=0.02,
                       radius=min(3, 10), rlen=nrow(indata)*10,
                       inv.alp.c=nrow(indata)*10/100)

    co.so.code <- co.so$code
  } else
  {
    co.so <- som.result
    co.so.code <- metadata
  }

  par(mfrow=c(co.dim, co.dim))
  par(mar=c(0, 0, 0, 0))

  for (j in (co.dim-1):0)
  {
    for (i in 0:(co.dim-1))
    {
      ij <- (i+1) + j * co.dim


      if (length(unique(group.labels)) > 1)
      {
        par(mar=c(5.4, 0, 0, 0))
        image(cbind(1:length(group.colors)), col=group.colors, axes=F)
        par(new=T, mar=c(0, 0, 0, 0))
      }

      plot(co.so.code[ij,], ylim=c(min(co.so.code),max(co.so.code)*1.2),
           xlab="", ylab="", axes=F, type="l", lwd=1.4)

      points(co.so.code[ij,], pch=16, cex=1)
      abline(h=0, col="darkgray", lty=2, lwd=1)

      par(new=T)
      plot(0, type="n", axes=F, xlab="", ylab="")
      text(0.65, 0.5, co.so$code.sum[ij, "nobs"], cex=1, adj=0)

      box()
    }
  }

  ### Plot SOM 10x10 grid, each group one line ###

  if (length(unique(group.labels)) > 1)
  {
    par(mfrow=c(co.dim, co.dim))
    par(mar=c(0, 0, 0, 0))

    for (j in (co.dim-1):0)
    {
      for (i in 0:(co.dim-1))
      {
        ij <- (i+1)+j*co.dim

        plot(co.so.code[ij,], ylim=c(min(co.so.code),max(co.so.code)*1.2),
             xlab="", ylab="", axes=F, type="n", lwd=1.4)

        par(new=T)
        plot(0, type="n", axes=F, xlab="", ylab="")
        text(0.65, 0.5, co.so$code.sum[ij, "nobs"], cex=1, adj=0)

        par(new=T)

        plot(0,
             xlim=c(1, max(table(group.labels))),
             ylim=c(min(co.so.code),max(co.so.code)*1.2),
             xlab="", ylab="", axes=F, type="n", lwd=1.4)

        for (i in 1:length(unique(group.labels)))
        {
          lines(co.so.code[ij, which(group.labels==unique(group.labels)[i])],
                col=groupwise.group.colors)[i], lwd=2)
        }

        abline(h=0, col="darkgray", lty=2, lwd=1)
        box()
      }
    }
  }

  dev.off()
}
