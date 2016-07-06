modules.chromosomes <- function(spot.list, main, path)
{
	if(is.null(chromosome.list))
  {
		return()
	}
	
  pdf(path, 21/2.54, 29.7/2.54)
 
  sorted.unique.gene.positions <- unlist( sapply( 1:length(chromosome.list), function(i) paste( names(chromosome.list)[i], names(chromosome.list[[i]]), sep=" " ) ) )
  sorted.unique.gene.positions <- sort.label( sorted.unique.gene.positions )
  n.all.chr.genes <- length(unlist(chromosome.list))

  intersect.counts <- matrix(0,
                             length(sorted.unique.gene.positions),
                             length(spot.list$spots),
                             dimnames=list(sorted.unique.gene.positions,
                                           names(spot.list$spots)))

  for (m in names(spot.list$spots))
  {
    intersect.counts.spot <- sapply(chromosome.list, function(x)
    {
      sapply(x, function(y)
      {
        length(intersect(spot.list$spots[[m]]$genes, y)) / length(y)
      })
    })

    intersect.counts.spot <- unlist(intersect.counts.spot)
    names(intersect.counts.spot) <- sub(".", " ", names(intersect.counts.spot), fixed=TRUE)
    intersect.counts[names(intersect.counts.spot), m] <- intersect.counts.spot
  }

  # write.csv2(intersect.counts, file=file.path(dirnames["csv"], paste(main, ".csv", sep="")))

  colkey <- c(colorRampPalette(c(rep("white",2),
                                 "green",
                                 rep("yellow2",1),
                                 rep("orange",2),
                                 rep("red",2),
                                 rep("darkred",2),
                                 rep("purple",1)))(2000),
              colorRampPalette(c("purple","black"))(3000),
              rep("black",5000))

  chr <- sapply(strsplit(rownames(intersect.counts), " ", fixed=TRUE), head, 1)
  chr.sep.level <- table(chr)[unique(chr)]
  chr.sep.level <- rev(chr.sep.level)
  chr.sep.level <- chr.sep.level[-length(chr.sep.level)]
  chr.sep.level <- sapply(seq_along(chr.sep.level), function(i) { sum(chr.sep.level[1:i]) })

  chr.lab.level <- table(chr)[unique(chr)]
  chr.lab.level <- rev(chr.lab.level)

  chr.lab.level <- c(0, sapply(1:(length(chr.lab.level)-1),
                               function(i) { sum(chr.lab.level[1:i]) })) + chr.lab.level/2

  par(mfrow=c(1,1), mar=c(2,3,2,2))

  image(1:ncol(intersect.counts),
        1:nrow(intersect.counts),
        t(intersect.counts[nrow(intersect.counts):1,]),
        axes=FALSE, xlab="", ylab="", col=colkey, main=paste("Chromosome Map,",main), zlim=c(0,1))

  title(main="(% band genes found in a spot)",line=0.1,cex.main=0.6)
  box()
  abline(h=chr.sep.level+0.5, lwd=0.1, lty=3)

  axis(1, 1:ncol(intersect.counts), labels = colnames(intersect.counts),
       las = 1, line = -0.5, tick = 0, cex.axis=0.8)

  axis(2, 1:nrow(intersect.counts), labels = rev(rownames(intersect.counts)),
       las = 2, line = -0.5, tick = 0, cex.axis=0.24)

  axis(2, chr.lab.level,  rev(unique(chr)), las=2, tick=FALSE, line=0.6)

  par(new=TRUE, mar=c(1,0,0,0))
  layout(matrix(c(0,0,0,0,1,0,0,0,0), 3, 3), c(1.6, 0.03, 0.005), c(0.115, 1.15, 2.15))
  image(matrix(1:1000, 1, 1000), col = colkey, axes=FALSE)

  axis(2, at=seq(0, 1, 0.1), paste(seq(0, 100, 10),"%",sep=""),
       las=2, tick=FALSE, pos=0.5, cex.axis=0.6)

  box()

  for (m in seq_along(spot.list$spots))
  {
    intersect.counts <- sapply(chromosome.list, function(x)
    {
      sapply(x, function(y)
      {
        length(intersect(spot.list$spots[[m]]$genes,  y)) / length(y)
      })
    })
    names(intersect.counts) <- names(chromosome.list)

    l <- c(rep(1,ceiling(length(chromosome.list)/2)), 2:(length(chromosome.list)+1))
    l <- c(l, rep(0, (ceiling(length(chromosome.list)/2)*3) - length(l)))
    layout(matrix(l, 3, ceiling(length(chromosome.list)/2), byrow=TRUE), heights=c(0.05,1,1))

    par(mar=c(0,0,0,0))
    plot(0, type="n", xlab="", ylab="", axes=FALSE, xlim=c(0,1))
    text(0.005, 0, paste("Spot",names(spot.list$spots)[m]), cex=2)
    text(0.2, 0, "(% band genes found in this spot)", cex=1.2)

    par(mar=c(1.5,2,0.2,0.3))
    for (chromosome in sort.label(names(chromosome.list)))
    {
      x <- intersect.counts[[chromosome]]

      image(1, seq_along(x),  matrix(rev(x), 1, length(x)), zlim=c(0, 1),
            axes=FALSE, xlab="", ylab="", col=colkey)

      box()
      axis(1, 1, labels = chromosome, las = 1, line = -0.8, tick = 0, cex.axis=1.2)
      axis(2, seq_along(x), labels = rev(names(x)), las = 2, line = -0.5, tick = 0, cex.axis=1)
    }
  }

  dev.off()
}
