pipeline.calcStatistics <- function()
{
  util.info("Calculating Single Gene Statistic")
 
  WAD.g.m <<- matrix(NA, nrow(indata), ncol(indata), dimnames=list(rownames(indata), colnames(indata)))
  
  for (m in 1:ncol(indata))
  {
    delta.e.g.m <- indata[,m]
  
    w.g.m <- (delta.e.g.m - min(delta.e.g.m)) / (max(delta.e.g.m) - min(delta.e.g.m))
    WAD.g.m[,m] <<- w.g.m * delta.e.g.m
  }

  
  # Calculate T-score and significance
  
  sd.g.m <- matrix(NA, nrow(indata), ncol(indata), dimnames=list(rownames(indata), colnames(indata)))
  
  t.g.m <<- matrix(NA, nrow(indata), ncol(indata), dimnames=list(rownames(indata), colnames(indata)))
  p.g.m <<- matrix(NA, nrow(indata), ncol(indata), dimnames=list(rownames(indata), colnames(indata)))
  
  n.0.m <<- rep(NA, ncol(indata))
  names(n.0.m) <<- colnames(indata)
  
  perc.DE.m <<- rep(NA, ncol(indata))
  names(perc.DE.m) <<- colnames(indata)
  
  fdr.g.m <<- matrix(NA, nrow(indata), ncol(indata), dimnames=list(rownames(indata), colnames(indata)))
  Fdr.g.m <<- matrix(NA, nrow(indata), ncol(indata), dimnames=list(rownames(indata), colnames(indata)))
  

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
  }

  ### Metagenes ###

  util.info("Calculating Metagene Statistic")

  t.m <<- p.m <<-
    matrix(NA, preferences$dim.1stLvlSom ^ 2, ncol(indata),
           dimnames=list(1:(preferences$dim.1stLvlSom ^ 2), colnames(indata)))

  t.m.help <- do.call(rbind, by(t.g.m, som.result$feature.BMU, colMeans))
  t.m[rownames(t.m.help),] <<- t.m.help

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
  }

}
