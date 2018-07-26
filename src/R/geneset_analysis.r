GeneSet.maxmean <- function(z, gs.def.list)
{
  s.plus <- pmax(z,0)
  s.minus <- pmin(z,0)
  
  S <- sapply(gs.def.list, function(x, S.plus, S.minus )
  {
    S.plus <- mean( s.plus[x$Genes] )
    S.minus <- mean( s.minus[x$Genes] )
    
    if( S.plus >= abs(S.minus) ) S.plus else S.minus
    
  }, s.plus, s.minus )
  
  return(S)
}


GeneSet.Fisher <- function(list.ids, all.ids, gs.def.list, sort=FALSE)
{
  list.ids <- list.ids[!is.na(list.ids)&!list.ids==""]

  fn <- function(x, list.ids, all.ids)
  {
    set.ids <- x$Genes

    n.list.set <- sum(set.ids %in% list.ids)
    n.set <- sum(! set.ids %in% list.ids)
    n.list <- sum(! list.ids %in% set.ids)
    n.none <- length(all.ids) - (n.list.set + n.list + n.set)

    return(fisher.test(matrix(c(n.list.set,n.set,n.list,n.none),2), alternative="greater")$p.value)
  }

  p.values <- sapply(gs.def.list, fn, list.ids, all.ids)
  names(p.values) <- names(gs.def.list)
  p.values[which(p.values > 1)] <- 1
  p.values[which(p.values < 1e-99)] <- 1e-99

  o <- c(seq_along(p.values))

  if (sort)
  {
    o <- order(p.values)
  }

  return(p.values[o])
}

GeneSet.GSZ <- function(list.symbols, all.gene.statistic, gs.def.list, sort=FALSE, N.min.list=10, N.min.set=10)
{
  N <- length(all.gene.statistic)
  N.list <- length(list.symbols)

  T.list <- sum(all.gene.statistic[list.symbols])
  mean.t <- mean(all.gene.statistic[list.symbols])
  var.t <- (1 / N.list) * sum((all.gene.statistic[list.symbols] - mean.t)^2)

  mean.t.all <- mean(all.gene.statistic)
  var.t.all <- (1 / length(all.gene.statistic)) * sum((all.gene.statistic - mean(all.gene.statistic))^2)

  all.gene.statistic <- all.gene.statistic[list.symbols]

  all.list.symbol.dummy <- rep(1,length(list.symbols))
  names(all.list.symbol.dummy) <- list.symbols

  GSZ <- sapply(gs.def.list, function(x, all.list.symbol.dummy,
                                      all.gene.statistic, N, N.list, T.list,
                                      mean.t, var.t, mean.t.all, var.t.all,
                                      N.min.list, N.min.set)
  {
    set.symbols <- x$Genes

    N.set <- length(set.symbols)
    N.plus <- sum(all.list.symbol.dummy[set.symbols], na.rm=TRUE)
    #exact: N.plus <- length(intersect(list.symbols, set.symbols))

    T.plus <- sum(all.gene.statistic[set.symbols], na.rm=TRUE)
    #exact: T.plus <- sum(all.gene.statistic[intersect(list.symbols, set.symbols)])

    T.minus <- sum(all.gene.statistic) - T.plus
    #exact: T.minus <- sum(all.gene.statistic[setdiff(list.symbols, set.symbols)])

    T.list <- T.plus + T.minus
    delta.T.list <- T.plus - T.minus

    E.N.plus <- N.list * N.set / N
    E.delta.T.list <- mean.t * (2 * E.N.plus - N.list)

    # SE exact:
    # SE.delta.T.list.2 <-
    #   4 * (var.t/(N.list-1) * (E.N.plus * (N.list - E.N.plus) - E.N.plus *
    #   (1 - N.set/N) * ((N-N.list)/(N-1))) + mean.t^2 * E.N.plus *
    #   (1 - N.set/N) * ((N-N.list)/(N-1)))

    # SE approximated:
    SE.delta.T.list.2 <-
      4 * E.N.plus * (1 - N.set/N) * (var.t + mean.t^2 * (1 - N.list/N))

    lambda <- 1 - min(1, sqrt(N.min.list / N.list * N.min.set / N.set))
    SE.0.2 <- 4 * E.N.plus * (1 - N.min.set/N) * (var.t.all + mean.t.all^2 * (1 - N.min.list/N))

    return((delta.T.list - E.delta.T.list) / sqrt(lambda * SE.delta.T.list.2 + (1 - lambda) * SE.0.2))

  }, all.list.symbol.dummy, all.gene.statistic, N, N.list, T.list, mean.t, var.t,
     mean.t.all, var.t.all, N.min.list, N.min.set)

  names(GSZ) <- names(gs.def.list)
  GSZ[which(is.infinite(GSZ))] <- NA

  if (sort)
  {
    o <- order(GSZ, decreasing=TRUE)
    GSZ <- GSZ[o]
  }

  return(GSZ)
}
