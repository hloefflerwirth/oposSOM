modules.CSV.sheets <- function(spot.list, main, path)
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
      gene <- indata[x,]
      return(suppressWarnings(cor(gene, spot.list$spotdata[m,])))
    })

    r.t <- r.genes / sqrt((1-r.genes^2) / (ncol(indata)-2))
    r.p <- 1 - pt(r.t, ncol(indata)-2)


    e.max <- apply(indata[spot.list$spots[[m]]$genes, ,drop=FALSE], 1, max)
    e.min <- apply(indata[spot.list$spots[[m]]$genes, ,drop=FALSE], 1, min)

    if (main %in% c("Sample-Underexpression"))
    {
      o <- names(sort(e.min, decreasing=FALSE))
    }  else
    {
      o <- names(sort(e.max, decreasing=TRUE))
    }

    out <- data.frame(Rank=c(seq_along(spot.list$spots[[m]]$genes)),
                      ID=o,
                      Symbol=gene.info$names[o])

    out <- cbind(out,
                 "mean expression"=indata.gene.mean[o],
                 "SD"=apply(indata[o, ,drop=FALSE], 1, sd),
                 "max delta e"=apply(indata[o, ,drop=FALSE], 1, max),
                 "min delta e"=apply(indata[o, ,drop=FALSE], 1, min),
                 "correlation"=r.genes[o],
                 "->t.score"=r.t[o],
                 "->p.value"=paste(r.p[o],"     ."),
                 "Metagene"=gene.info$coordinates[o],
                 "Chromosome"=paste( gene.info$chr.name[o], gene.info$chr.band[o]),
                 "Description"=gene.info$descriptions[o])

    write.csv2(out, file.path(path, basename))
  }
}