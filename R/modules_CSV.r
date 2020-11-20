modules.CSV.sheets <- function(env, spot.list, main, path)
{
  for (m in seq_along(spot.list$spots))
  {
    basename <- paste(main, " ", names(spot.list$spots)[m], ".csv", sep="")

    if (length(spot.list$spots[[m]]$genes) <= 0 || length(spot.list$spots[[m]]$genes) > 2000)
    {
      next
    }

    ## CSV Table
    r.genes <- sapply(spot.list$spots[[m]]$genes, function(x)
    {
      gene <- env$indata[x,]
      return(suppressWarnings(cor(gene, spot.list$spotdata[m,])))
    })

    r.t <- r.genes / sqrt((1-r.genes^2) / (ncol(env$indata)-2))
    r.p <- 1 - pt(r.t, ncol(env$indata)-2)


    e.max <- apply(env$indata[spot.list$spots[[m]]$genes, ,drop=FALSE], 1, max)
    e.min <- apply(env$indata[spot.list$spots[[m]]$genes, ,drop=FALSE], 1, min)

    if (main %in% c("Underexpression Spots"))
    {
      o <- names(sort(e.min, decreasing=FALSE))
    }  else
    {
      o <- names(sort(e.max, decreasing=TRUE))
    }

    out <- data.frame(Rank=c(seq_along(spot.list$spots[[m]]$genes)),
                      ID=o,
                      Symbol=env$gene.info$names[o])

    out <- cbind(out,
                 "mean expression"=env$indata.gene.mean[o],
                 "SD"=apply(env$indata[o, ,drop=FALSE], 1, sd),
                 "max delta e"=apply(env$indata[o, ,drop=FALSE], 1, max),
                 "min delta e"=apply(env$indata[o, ,drop=FALSE], 1, min),
                 "correlation"=r.genes[o],
                 "->t.score"=r.t[o],
                 "->p.value"=paste(r.p[o],"     ."),
                 "Metagene"=env$gene.info$coordinates[o],
                 "Chromosome"=paste( env$gene.info$chr.name[o], env$gene.info$chr.band[o]),
                 "Description"=env$gene.info$descriptions[o])

    env$csv.function(out, file.path(path, basename))
  }
}