pipeline.groupAssignment <- function()
{
  if (length(unique(group.labels)) < 2)
  {
    group.bootstrap.score <<- rep(100, ncol(indata))
    return()
  }

  group.assignment <- match(group.labels, unique(group.labels))
  names(group.assignment) <- colnames(indata)

  bootstrap.error <- rep(0,ncol(indata))
  names(bootstrap.error) <- colnames(indata)

  n.bootstrap <- min(1000, 10*ncol(indata))

  # Note:
  #
  # This parallelization tends to be (a lot) slower than non-parallelized code,
  # when working on a small number of samples (ncol(indata) < 40 ?), but
  # gives great advantages for big datasets.
  cl <- makeCluster(preferences$max.parallel.cores)

  bootstrap.error <- parSapply(cl, 1:n.bootstrap, function(i)
  {
    resample <- sample(colnames(indata), ncol(indata), replace=T)

    suppressWarnings({
      km <- kmeans(t(metadata[,resample]),
                   centers=t(group.metadata),
                   iter.max=10,
                   algorithm = "Lloyd")
    })

    if (length(km) > 0)
    {
      return(unique(names(which(km$cluster[resample] != group.assignment[resample]))))
    }
  })

  stopCluster(cl)

  bootstrap.error <- table(unlist(bootstrap.error))[colnames(indata)]
  names(bootstrap.error) <- colnames(indata)
  bootstrap.error[which(is.na(bootstrap.error))] <- 0

  group.bootstrap.score <<- round(100*(1 - bootstrap.error / n.bootstrap), 1)
}
