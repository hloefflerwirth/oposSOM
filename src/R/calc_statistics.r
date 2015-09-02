pipeline.calcStatistics <- function()
{
  verbose = (output.paths["LPE"] != "")


  
  util.info("Calculating Single Gene Statistic")
  util.progress(0, 48)
 
  WAD.g.m <<- matrix(NA, nrow(indata), ncol(indata), dimnames=list(rownames(indata), colnames(indata)))
  
  for (m in 1:ncol(indata))
  {
    delta.e.g.m <- indata[,m]
  
    w.g.m <- (delta.e.g.m - min(delta.e.g.m)) / (max(delta.e.g.m) - min(delta.e.g.m))
    WAD.g.m[,m] <<- w.g.m * delta.e.g.m
  }

  
  # Calculate T-score and significance
  progress.max <- ncol(indata)
  progress.current <- 0

  
  sd.g.m <- matrix(NA, nrow(indata), ncol(indata), dimnames=list(rownames(indata), colnames(indata)))
  
  t.g.m <<- matrix(NA, nrow(indata), ncol(indata), dimnames=list(rownames(indata), colnames(indata)))
  p.g.m <<- matrix(NA, nrow(indata), ncol(indata), dimnames=list(rownames(indata), colnames(indata)))
  
  n.0.m <<- rep(NA, ncol(indata))
  names(n.0.m) <<- colnames(indata)
  
  perc.DE.m <<- rep(NA, ncol(indata))
  names(perc.DE.m) <<- colnames(indata)
  
  fdr.g.m <<- matrix(NA, nrow(indata), ncol(indata), dimnames=list(rownames(indata), colnames(indata)))
  Fdr.g.m <<- matrix(NA, nrow(indata), ncol(indata), dimnames=list(rownames(indata), colnames(indata)))
  

  {
    o <- order(indata.gene.mean)
    sdo <- apply(indata, 1, sd)[o]
    col <- Get.Running.Average(sdo, min(200, round(nrow(indata) * 0.02)))
    col[which(is.nan(col))] <- 0.0000000001
    col[which(col == 0)] <- 0.0000000001

    for (i in seq(length(col)-1, 1))
    {
      col[i] <- max(col[i], col[i+1])
    }

    sd.g.m[o,] <- col

    t.g.m <<- apply(indata, 2, function(x, root)
    {
      return(root * x / sd.g.m[,1])
    }, sqrt(ncol(indata)))

    progress.current <- progress.current + (0.6 * ncol(indata))
    util.progress(progress.current, progress.max)
  } 



  ### calculate significance and fdr ###

  for (m in 1:ncol(indata))
  {
#    p.g.m[,m] <<- 2 - 2*pt( abs(t.g.m[,m]), ncol(indata) - 1 )
    
    suppressWarnings({
      try.res <- try({
#        fdrtool.result <- fdrtool(p.g.m[,m], statistic="pvalue", verbose=FALSE, plot=FALSE)
        fdrtool.result <- fdrtool(t.g.m[,m], verbose=FALSE, plot=FALSE)
      }, silent=TRUE)
    })

    if (class(try.res) != "try-error")
    {
      p.g.m[,m] <<- fdrtool.result$pval
      fdr.g.m[,m] <<- fdrtool.result$lfdr
      Fdr.g.m[,m] <<- fdrtool.result$qval

      n.0.m[m] <<- fdrtool.result$param[1,"eta0"]
      perc.DE.m[m] <<- 1 - n.0.m[m]
    } else # happens for eg phenotype data
    {
      p.g.m[,m] <<- order(indata[,m]) / nrow(indata)
      fdr.g.m[,m] <<- p.g.m[,m]
      Fdr.g.m[,m] <<- p.g.m[,m]

      n.0.m[m] <<- 0.5
      perc.DE.m[m] <<- 1 - n.0.m[m]
    }

    progress.current <- progress.current + 0.4
    util.progress(progress.current, progress.max)
  }

  util.progress.terminate()


  ### error plots ###

  if (verbose)
  {
    filename <- file.path(output.paths["LPE"], "all_sample_LPE.bmp")
    util.info("Writing:", filename)

    bmp(filename, 600, 600)
    par(mar=c(5, 6, 4, 5))

    plot(apply(indata, 1, sd) ~ indata.gene.mean,
         xlab=expression(e[g]),
         ylab="",
         main="Locally pooled error estimate (LPE)",
         las=1,
         cex.main=1.5,
         cex.lab=2,
         cex.axis=2)      
    
    mtext(expression(sigma[g]), side=2, line=4, las=2, cex=2)
    points(sd.g.m[,1] ~ indata.gene.mean, col="green", pch=16)
    legend("topright","LPE",lwd=4,col="green")
    dev.off()
  }


  ### Metagenes ###

  progress.current <- 0
  progress.max <- ncol(indata)

  util.info("Calculating Metagene Statistic")
  util.progress(progress.current, progress.max)

  t.m <<- p.m <<-
    matrix(NA, preferences$dim.1stLvlSom ^ 2, ncol(indata),
           dimnames=list(1:(preferences$dim.1stLvlSom ^ 2), colnames(indata)))

  t.m.help <- do.call(rbind, by(t.g.m, som.nodes, colMeans))
  t.m[rownames(t.m.help),] <<- t.m.help

  progress.current <- 0.4 * progress.max
  util.progress(progress.current, progress.max)

  for (m in 1:ncol(indata))
  {
    suppressWarnings({
      try.res <- try({
        fdrtool.result <- fdrtool(as.vector(na.omit(t.m[,m])), verbose=FALSE, plot=FALSE)
      }, silent=TRUE)
    })

    if (class(try.res) != "try-error")
    {
      p.m[which(!is.na(t.m[,m])),m] <<- fdrtool.result$pval
    } else # happens for eg phenotype data
    {
      p.m[which(!is.na(t.m[,m])),m] <<- t.m[which(!is.na(t.m[,m])),m] / max(t.m[,m], na.rm=TRUE)
    }

    progress.current <- progress.current + 0.6
    util.progress(progress.current, progress.max)
  }

  util.progress.terminate()

}
