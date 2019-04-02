Get.Running.Average = function(v, n=length(v)/100)
{

  ret = c(NA, length(v))


  window = sum (v[1 : (2 * n + 1)])

  for (i in seq_along(v))
  {
    if (i > n + 1 &&
        i <= length(v) - n)
                             {
      window = window + v[i + n]
      window = window - v[i - n - 1]
    }

    ret[i] = window / (2 * n + 1)

  }

  names(ret) = names(v)
  return(ret)

}

Smooth.Matrix <- function (m, window = length(v)/100) 
{
  bins <- ceiling(1+window/2):floor(nrow(m)-window/2)
  
  t( sapply( bins, function(i)
  {
    colMeans( m[ c( (i-window/2) : (i+window/2 )), ] )
  } ) )
}
