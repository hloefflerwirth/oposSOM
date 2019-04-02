Quantile.Normalization = function(data)
{
  quantiles = c(rep(0, nrow(data)))
  na.cols = c(rep(0, nrow(data)))

  for (i in 1:ncol(data))
  {  
    data.col = data [order(data [,i], decreasing=TRUE) ,i]
    not.na = which(!is.na(data.col))
    
    quantiles[not.na] = quantiles[not.na] + data.col[not.na]
    na.cols[not.na] = na.cols[not.na] + 1
  }
  quantiles = quantiles / na.cols
  
  for (i in 1:ncol(data))
  {
    data [order(data [,i], decreasing=TRUE) ,i] = quantiles 
  }

  return(invisible (data))
}
