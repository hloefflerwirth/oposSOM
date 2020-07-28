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

Sample.GSZ <- function(gene.set,t.ensID.m,mean.t.all,sd.t.all)
{
  mean.t.set <- colMeans( t.ensID.m[gene.set$Genes,] )
  
  GSZ <- sqrt(length(gene.set$Genes)) * ( mean.t.set - mean.t.all ) / sd.t.all
}