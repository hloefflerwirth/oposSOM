pipeline.diffExpressionStatistics <- function(env)
{
  util.info("Calculating Single Gene Statistic")
  
  ### Single Genes ###
  
  env$p.g.m <- chunk.apply.rows( env$indata, function(x)
  {
    sapply( seq(x), function(m)
    {
      if( all(x[-m] == x[-m][1]) ) return(1)
      
      return( t.test( x[m], x[-m], var.equal=TRUE )$p.value )
    })
  } )
  colnames(env$p.g.m) <- colnames(env$indata)
  
  
  env$n.0.m <- rep(NA, ncol(env$indata))
  names(env$n.0.m) <- colnames(env$indata)
  
  env$perc.DE.m <- rep(NA, ncol(env$indata))
  names(env$perc.DE.m) <- colnames(env$indata)
  
  env$fdr.g.m <- matrix(NA, nrow(env$indata), ncol(env$indata), dimnames=list(rownames(env$indata), colnames(env$indata)))
  

  for (m in 1:ncol(env$indata))
  {
    suppressWarnings({ try.res <- try({
       fdrtool.result <- fdrtool(env$p.g.m[,m], statistic="pvalue", verbose=FALSE, plot=FALSE)
    }, silent=TRUE) })

    if (!is(try.res,"try-error"))
    {
      env$fdr.g.m[,m] <- fdrtool.result$lfdr

      env$n.0.m[m] <- fdrtool.result$param[1,"eta0"]
      env$perc.DE.m[m] <- 1 - env$n.0.m[m]
    } else # happens for eg phenotype data
    {
      env$p.g.m[,m] <- order(env$indata[,m]) / nrow(env$indata)
      env$fdr.g.m[,m] <- env$p.g.m[,m]

      env$n.0.m[m] <- 1
      env$perc.DE.m[m] <- 0
    }

  }

  ### Metagenes ###

  util.info("Calculating Metagene Statistic")

  env$p.m <- matrix(NA, env$preferences$dim.1stLvlSom ^ 2, ncol(env$indata),
           dimnames=list(1:(env$preferences$dim.1stLvlSom ^ 2), colnames(env$indata)))

  p.m.help <- do.call(rbind, by(env$p.g.m, env$som.result$feature.BMU, colMeans))
  env$p.m[rownames(p.m.help),] <- p.m.help
  
  return(env)
}
