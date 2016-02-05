pipeline.2ndLvlModuleCorrelation <- function(s, hcl)
{
  pcm <- cor(s)
  pcm <- pcm[hcl$order,hcl$order]
  core.rect.thresh <- quantile(pcm, 0.9)
  pcm.mask <- matrix(NA,ncol(indata),ncol(indata),dimnames=dimnames(pcm))

  for (core.rect.dim in ncol(indata)*c(0.1,0.05,0.02))
  {
    x.i <- 1
    y.i <- 1

    core.rect.sum.thresh <- core.rect.thresh * core.rect.dim^2

    for (x.i in 1:(ncol(indata)-core.rect.dim+1))
    {
      for (y.i in x.i:(ncol(indata)-core.rect.dim+1))
      {
        core.sum <- sum(pcm[c(x.i:(x.i+core.rect.dim-1)), c(y.i:(y.i+core.rect.dim-1))])

        if (core.sum > core.rect.sum.thresh)
        {
          pcm.mask[x.i:(x.i+core.rect.dim-1), y.i:(y.i+core.rect.dim-1)] <- 1
          pcm.mask[y.i:(y.i+core.rect.dim-1), x.i:(x.i+core.rect.dim-1)] <- 1
        }  else if (core.sum < -core.rect.sum.thresh)
        {
          pcm.mask[x.i:(x.i+core.rect.dim-1), y.i:(y.i+core.rect.dim-1)] <- -1
          pcm.mask[y.i:(y.i+core.rect.dim-1), x.i:(x.i+core.rect.dim-1)] <- -1
        }
      }
    }
  }

  pcm.mask[which(is.na(pcm.mask))] <- 0

  if (any(pcm.mask != 0))
  {
    boxes.pos <- list()

    while (nrow(which(pcm.mask==1,arr.ind=TRUE)) > 0)
    {
      x.i <- which(pcm.mask==1,arr.ind=TRUE)[1,1]
      y.i <- which(pcm.mask==1,arr.ind=TRUE)[1,2]

      dim.x <- 1
      dim.y <- 1

      while ((x.i+dim.x) <= ncol(indata) &&
             (y.i+dim.y) <= ncol(indata) &&
             all(pcm.mask[x.i:(x.i+dim.x), y.i:(y.i+dim.y)] == 1, na.rm=TRUE))
      {
        dim.x <- dim.x + 1
        dim.y <- dim.y + 1
      }

      grown <- TRUE

      while (grown)
      {
        grown <- FALSE

        if ((x.i+dim.x) <= ncol(indata) &&
            all(pcm.mask[x.i:(x.i+dim.x), y.i:(y.i+dim.y-1)] == 1, na.rm=TRUE))
        {
          dim.x <- dim.x + 1
          grown <- TRUE
        } else if ((y.i+dim.y) <= ncol(indata) &&
                   all(pcm.mask[x.i:(x.i+dim.x-1), y.i:(y.i+dim.y)] == 1, na.rm=TRUE))
        {
          dim.y <- dim.y + 1
          grown <- TRUE
        }
      }

      boxes.pos[[length(boxes.pos)+1]] <- c(x=x.i, y=y.i, dim.x=dim.x, dim.y=dim.y)
      pcm.mask[x.i:(x.i+dim.x-1), y.i:(y.i+dim.y-1)] <- 0
    }

    boxes.neg <- list()

    while (nrow(which(pcm.mask==-1,arr.ind=TRUE)) > 0)
    {
      x.i <- which(pcm.mask==-1,arr.ind=TRUE)[1,1]
      y.i <- which(pcm.mask==-1,arr.ind=TRUE)[1,2]

      dim.x <- 1
      dim.y <- 1

      while ((x.i+dim.x) <= ncol(indata) &&
             (y.i+dim.y) <= ncol(indata) &&
             all(pcm.mask[x.i:(x.i+dim.x), y.i:(y.i+dim.y)] == -1, na.rm=TRUE))
      {
        dim.x <- dim.x + 1
        dim.y <- dim.y + 1
      }

      grown <- TRUE

      while (grown)
      {
        grown <- FALSE

        if ((x.i+dim.x) <= ncol(indata) &&
            all(pcm.mask[x.i:(x.i+dim.x), y.i:(y.i+dim.y-1)] == -1, na.rm=TRUE))
        {
          dim.x <- dim.x + 1
          grown <- TRUE
        } else if ((y.i+dim.y) <= ncol(indata) &&
                   all(pcm.mask[x.i:(x.i+dim.x-1), y.i:(y.i+dim.y)] == -1, na.rm=TRUE))
        {
          dim.y <- dim.y + 1
          grown <- TRUE
        }
      }

      boxes.neg[[length(boxes.neg)+1]] <- c(x=x.i, y=y.i, dim.x=dim.x, dim.y=dim.y)
      pcm.mask[x.i:(x.i+dim.x-1), y.i:(y.i+dim.y-1)] <- 0
    }

    boxes.pos <- boxes.pos[which(sapply(boxes.pos, function(x) { x['dim.x'] * x['dim.y'] }) >
                                 ceiling((ncol(indata)*0.02)^2))]

    boxes.neg <- boxes.neg[which(sapply(boxes.neg, function(x) { x['dim.x'] * x['dim.y'] }) >
                                 ceiling((ncol(indata)*0.02)^2))]

    boxes.pos <- boxes.pos[which(sapply(boxes.pos, function(x) { min(x['dim.x'], x['dim.y']) }) >
                                 ceiling(ncol(indata)*0.01))]

    boxes.neg <- boxes.neg[which(sapply(boxes.neg, function(x) { min(x['dim.x'], x['dim.y']) }) >
                                 ceiling(ncol(indata)*0.01))]




    #### Plot PCM Clone with Spot Correlation Modules ####
    layout(matrix(c(0, 0, 4, 0, 0, 1, 5, 2, 3), 3),
           widths=c(1, 0.2, 4), heights=c(1.2, 0.2, 4), respect=TRUE)

    par(mar = c(8, 0, 0, 0.5))
    image(rbind(1:ncol(indata)), col=group.colors[hcl$order], axes=FALSE)

    par(mar = c(0.5, 0, 0, 8))
    image(cbind(1:ncol(indata)), col=group.colors[hcl$order], axes=FALSE)

    par(mar = c(8, 0, 0, 8))

    image(x=c(0:ncol(indata)), y=c(0:ncol(indata)), z=pcm,col=color.palette.heatmaps(1000),
          axes=FALSE, xlab="", ylab="")

    box()

    if (length(boxes.pos) > 0)
    {
      dummy <- sapply(boxes.pos, function(x)
      {
        rect(x['x']-1, x['y']-1, x['x']+x['dim.x']-1, x['y']+x['dim.y']-1, col="red", border="black")
      })
    }

    if (length(boxes.neg) > 0)
    {
      dummy <- sapply(boxes.neg, function(x)
      {
        rect(x['x']-1, x['y']-1, x['x']+x['dim.x']-1, x['y']+x['dim.y']-1, col="blue", border="black")
      })
    }

    spotdata.sd.scale <- 1

    if (length(boxes.pos) > 0)
    {
      dummy <- sapply(boxes.pos, function(x)
      {
        samples.x <- colnames(pcm)[c(x['x']:(x['x']+x['dim.x']-1))]
        samples.y <- colnames(pcm)[c(x['y']:(x['y']+x['dim.y']-1))]

        spots.x <- names(which(rowMeans(spot.list.dmap$spotdata[,samples.x, drop=FALSE]) >
                               sd(as.vector(spot.list.dmap$spotdata)) * spotdata.sd.scale))

        if (length(spots.x)==0)
        {
          spots.x <- names(which.max(rowMeans(spot.list.dmap$spotdata[,samples.x, drop=FALSE])))
        }

        spots.y <- names(which(rowMeans(spot.list.dmap$spotdata[,samples.y, drop=FALSE]) >
                               sd(as.vector(spot.list.dmap$spotdata)) * spotdata.sd.scale))

        if (length(spots.y)==0)
        {
          spots.y <- names(which.max(rowMeans(spot.list.dmap$spotdata[,samples.y, drop=FALSE])))
        }

        text(mean(c(x['x']:(x['x']+x['dim.x'])))-1, mean(c(x['y']:(x['y']+x['dim.y'])))-1,
             paste(intersect(spots.x, spots.y), collapse=","), cex=0.8, col="gray90")
      })
    }

    if (length(boxes.neg) > 0)
    {
      dummy <- sapply(boxes.neg, function(x)
      {
        samples.x <- colnames(pcm)[c(x['x']:(x['x']+x['dim.x']-1))]
        samples.y <- colnames(pcm)[c(x['y']:(x['y']+x['dim.y']-1))]

        spots.x1 <- names(which(rowMeans(spot.list.dmap$spotdata[,samples.x, drop=FALSE]) >
                                sd(as.vector(spot.list.dmap$spotdata)) * spotdata.sd.scale))

        if (length(spots.x1) == 0)
        {
          spots.x1 <- names(which.max(rowMeans(spot.list.dmap$spotdata[,samples.x, drop=FALSE])))
        }

        spots.y1 <- names(which(rowMeans(spot.list.dmap$spotdata[,samples.y, drop=FALSE]) <
                                sd(as.vector(spot.list.dmap$spotdata)) * spotdata.sd.scale))

        if (length(spots.y1) == 0)
        {
          spots.y1 <- names(which.min(rowMeans(spot.list.dmap$spotdata[,samples.y, drop=FALSE])))
        }

        spots.x2 <- names(which(rowMeans(spot.list.dmap$spotdata[,samples.x, drop=FALSE]) <
                                sd(as.vector(spot.list.dmap$spotdata)) * spotdata.sd.scale))

        if (length(spots.x2) == 0)
        {
          spots.x2 <- names(which.min(rowMeans(spot.list.dmap$spotdata[,samples.x, drop=FALSE])))
        }

        spots.y2 <- names(which(rowMeans(spot.list.dmap$spotdata[,samples.y, drop=FALSE]) >
                                sd(as.vector(spot.list.dmap$spotdata)) * spotdata.sd.scale))

        if (length(spots.y2) == 0)
        {
          spots.y2 <- names(which.max(rowMeans(spot.list.dmap$spotdata[,samples.y, drop=FALSE])))
        }

        text(mean(c(x['x']:(x['x']+x['dim.x'])))-1, mean(c(x['y']:(x['y']+x['dim.y'])))-1,
             paste(c(intersect(spots.x1, spots.y1), intersect(spots.x2, spots.y2)), collapse=","),
             cex=0.8, col="gray90")
      })
    }

    ddr <- as.dendrogram(hcl)

    par(mar=c(8, 0, 0, 0))
    plot(ddr, horiz=TRUE, axes=FALSE, yaxs="i", leaflab="none")
    par(mar=c(0, 0, 1, 8))
    plot(ddr, axes=FALSE, xaxs="i", leaflab="none")
    title(paste("Expression module correlation map (D-clusters),",metagene.filter.list[[i]]$n), cex.main = 1.5)
  }
}
