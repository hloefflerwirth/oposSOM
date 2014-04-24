col.pix <- function(m, x, y, col, dim, threshold=NA)
{
  hash.list <- c(x + y * dim)
  pos.list <- list(c(x, y))

  while (length(hash.list) > 0)
  {
    x <- pos.list[[1]][1]
    y <- pos.list[[1]][2]

    if ((is.na(threshold) && !is.na(m[x,y]) && m[x,y] != col) ||
        (!is.na(threshold) && m[x,y] > threshold))
    {
      m[x,y] <- col

      if (x-1 > 0 && length(intersect((x-1) + (y) * dim, hash.list)) == 0)
      {
        hash.list <- c(hash.list, (x-1) + (y) * dim)
        pos.list[[length(pos.list)+1]] <- c(x-1, y)
      }

      if (y-1 > 0 && length(intersect((x) + (y-1) * dim, hash.list)) == 0)
      {
        hash.list <- c(hash.list, (x) + (y-1) * dim)
        pos.list[[length(pos.list)+1]] <- c(x, y-1)
      }

      if (x+1 <= nrow(m) && length(intersect((x+1) + (y) * dim, hash.list)) == 0)
      {
        hash.list <- c(hash.list, (x+1) + (y) * dim)
        pos.list[[length(pos.list)+1]] <- c(x+1, y)
      }

      if (y+1 <= ncol(m) && length(intersect((x) + (y+1) * dim, hash.list)) == 0)
      {
        hash.list <- c(hash.list, (x) + (y+1) * dim)
        pos.list[[length(pos.list)+1]] <- c(x, y+1)
      }
    }

    hash.list <- hash.list[-1]
    pos.list <- pos.list[-1]
  }

  return(m)
}

get.neighbors <- function(x, y, dim)
{
  ret <- list()

  if (!missing("y"))
  {
    if (x > 1)
    {
      ret[['l']] <- c(x-1, y)
    }

    if (x < dim)
    {
      ret[['r']] <- c(x+1, y)
    }

    if (y > 1)
    {
      ret[['u']] <- c(x, y-1)
    }

    if (y < dim)
    {
      ret[['d']] <- c(x, y+1)
    }
  }else
  {
    y <- (x-1) %/% dim+1
    x <- (x-1) %% dim+1

    if (x > 1)
    {
      ret[['l']] <- (y-1) * dim + x - 1
    }

    if (x < dim)
    {
      ret[['r']] <- (y-1) * dim + x + 1
    }

    if (y > 1)
    {
      ret[['u']] <- (y-2) * dim + x
    }

    if (y < dim)
    {
      ret[['d']] <- (y) * dim + x
    }
  }

  return(ret)
}

pipeline.detectSpotsSamples <- function()
{
  GS.infos.samples <<- list()

  for (j in 1:ncol(indata))
  {
    GS.infos.samples[[j]] <<- list()

    mask <- matrix(NA, preferences$dim.som1, preferences$dim.som1)
    blob <- matrix(metadata[,j], preferences$dim.som1, preferences$dim.som1)

    if (diff(sign(range(metadata[,j]))) != 0) # values lie in + and - regions
    {
      mask[which(blob > max(blob) * preferences$sample.spot.cutoff)] <- -1
      mask[which(blob < min(blob) * preferences$sample.spot.cutoff)] <- -2
    } else
    {
      mask[which(blob > quantile(blob,0.9))] <- -1
      mask[which(blob < quantile(blob,0.1))] <- -2
    }

    GS.infos.samples[[j]]$regulated <<- mask

    spot.i <- 0
    spot.updown <- c()

    while (nrow(which(mask == -1, arr.ind=T)) > 0)
    {
      start.pix <- which(mask == -1, arr.ind=T)[1,]
      spot.i <- spot.i + 1
      mask <- col.pix(mask, start.pix[1], start.pix[2], spot.i, preferences$dim.som1)
      spot.updown  <- c(spot.updown, "overexpressed")
    }

    while (nrow(which(mask == -2, arr.ind=T)) > 0)
    {
      start.pix <- which(mask == -2, arr.ind=T)[1,]
      spot.i <- spot.i + 1
      mask <- col.pix(mask, start.pix[1], start.pix[2], spot.i, preferences$dim.som1)
      spot.updown <- c(spot.updown, "underexpressed")
    }

    GS.infos.samples[[j]]$spots <<- list()

    if (spot.i > 0)
    {
      for (spot.ii in 1:spot.i)
      {
        GS.infos.samples[[j]]$spots[[spot.ii]] <<- list()
        GS.infos.samples[[j]]$spots[[spot.ii]]$type <<- spot.updown[spot.ii]
        GS.infos.samples[[j]]$spots[[spot.ii]]$metagenes <<- which(mask == spot.ii)

        GS.infos.samples[[j]]$spots[[spot.ii]]$genes <<-
          names(som.nodes)[which(som.nodes %in% which(mask == spot.ii))]

        GS.infos.samples[[j]]$spots[[spot.ii]]$mask <<-
          matrix(NA, preferences$dim.som1, preferences$dim.som1)

        GS.infos.samples[[j]]$spots[[spot.ii]]$mask[which(mask == spot.ii)] <<- 1
      }
    }
  }

  names(GS.infos.samples) <<- colnames(indata)
}
