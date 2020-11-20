pipeline.calcStatistics <- function(env)
{
  util.info("Calculating Single Gene Statistic")
 
  env$WAD.g.m <- matrix(NA, nrow(env$indata), ncol(env$indata), dimnames=list(rownames(env$indata), colnames(env$indata)))
  
  for (m in 1:ncol(env$indata))
  {
    delta.e.g.m <- env$indata[,m]
  
    w.g.m <- (delta.e.g.m - min(delta.e.g.m)) / (max(delta.e.g.m) - min(delta.e.g.m))
    env$WAD.g.m[,m] <- w.g.m * delta.e.g.m
  }

  
  # Calculate T-score and significance
  
  sd.g.m <- matrix(NA, nrow(env$indata), ncol(env$indata), dimnames=list(rownames(env$indata), colnames(env$indata)))
  
  env$t.g.m <- matrix(NA, nrow(env$indata), ncol(env$indata), dimnames=list(rownames(env$indata), colnames(env$indata)))
  env$p.g.m <- matrix(NA, nrow(env$indata), ncol(env$indata), dimnames=list(rownames(env$indata), colnames(env$indata)))
  
  env$n.0.m <- rep(NA, ncol(env$indata))
  names(env$n.0.m) <- colnames(env$indata)
  
  env$perc.DE.m <- rep(NA, ncol(env$indata))
  names(env$perc.DE.m) <- colnames(env$indata)
  
  env$fdr.g.m <- matrix(NA, nrow(env$indata), ncol(env$indata), dimnames=list(rownames(env$indata), colnames(env$indata)))
  env$Fdr.g.m <- matrix(NA, nrow(env$indata), ncol(env$indata), dimnames=list(rownames(env$indata), colnames(env$indata)))
  

  o <- order(env$indata.gene.mean)
  sdo <- apply(env$indata, 1, sd)[o]
  col <- Get.Running.Average(sdo, min(200, round(nrow(env$indata) * 0.02)))
  col[which(is.nan(col))] <- 0.0000000001
  col[which(col == 0)] <- 0.0000000001

  for (i in seq(length(col)-1, 1))
  {
    col[i] <- max(col[i], col[i+1])
  }

  sd.g.m[o,] <- col

  env$t.g.m <- apply(env$indata, 2, function(x, root)
  {
    return(root * x / sd.g.m[,1])
  }, sqrt(ncol(env$indata)))



  ### calculate significance and fdr ###

  for (m in 1:ncol(env$indata))
  {
#    p.g.m[,m] <<- 2 - 2*pt( abs(t.g.m[,m]), ncol(indata) - 1 )
    
    suppressWarnings({
      try.res <- try({
#        fdrtool.result <- fdrtool(p.g.m[,m], statistic="pvalue", verbose=FALSE, plot=FALSE)
        fdrtool.result <- fdrtool(env$t.g.m[,m], verbose=FALSE, plot=FALSE)
      }, silent=TRUE)
    })

    if (!is(try.res,"try-error"))
    {
      env$p.g.m[,m] <- fdrtool.result$pval
      env$fdr.g.m[,m] <- fdrtool.result$lfdr
      env$Fdr.g.m[,m] <- fdrtool.result$qval

      env$n.0.m[m] <- fdrtool.result$param[1,"eta0"]
      env$perc.DE.m[m] <- 1 - env$n.0.m[m]
    } else # happens for eg phenotype data
    {
      env$p.g.m[,m] <- order(env$indata[,m]) / nrow(env$indata)
      env$fdr.g.m[,m] <- env$p.g.m[,m]
      env$Fdr.g.m[,m] <- env$p.g.m[,m]

      env$n.0.m[m] <- 0.5
      env$perc.DE.m[m] <- 1 - env$n.0.m[m]
    }
  }

  ### Metagenes ###

  util.info("Calculating Metagene Statistic")

  env$t.m <- env$p.m <-
    matrix(NA, env$preferences$dim.1stLvlSom ^ 2, ncol(env$indata),
           dimnames=list(1:(env$preferences$dim.1stLvlSom ^ 2), colnames(env$indata)))

  t.m.help <- do.call(rbind, by(env$t.g.m, env$som.result$feature.BMU, colMeans))
  env$t.m[rownames(t.m.help),] <- t.m.help

  for (m in 1:ncol(env$indata))
  {
    suppressWarnings({
      try.res <- try({
        fdrtool.result <- fdrtool(as.vector(na.omit(env$t.m[,m])), verbose=FALSE, plot=FALSE)
      }, silent=TRUE)
    })

		if( !is(try.res,"try-error") )
    {
		  env$p.m[which(!is.na(env$t.m[,m])),m] <- fdrtool.result$pval
    } else # happens for eg phenotype data
    {
      env$p.m[which(!is.na(env$t.m[,m])),m] <- env$t.m[which(!is.na(env$t.m[,m])),m] / max(env$t.m[,m], na.rm=TRUE)
    }
  }
  return(env)
}
