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
      ret[['d']] <- c(x, y-1)
    }

    if (y < dim)
    {
      ret[['u']] <- c(x, y+1)
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
      ret[['d']] <- (y-2) * dim + x
    }

    if (y < dim)
    {
      ret[['u']] <- (y) * dim + x
    }
  }

  return(ret)
}
