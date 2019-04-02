get.beta.statistic = function(set.data, ref.vector, weights)
{
  if (missing("ref.vector"))
  {
    ref.vector = colMeans(set.data)
  }
  
  
  ref.mean = mean(ref.vector)
  
  
  
  sum.r = apply(set.data, 1, function(x)
  {
    #  t(ref.vector - ref.mean) %*% (x - mean(x)) / sqrt(t(ref.vector - ref.mean) %*% (ref.vector - ref.mean))      
    
    sum((ref.vector - ref.mean) * (x - mean(x))) / sqrt(sum((ref.vector - ref.mean) ^ 2))          
    
  })
  
  
  if (missing("weights"))
  {
    sum.r = sum(sum.r)
    
  } else
  {
    sum.r = 1/sum(weights) * sum(weights * sum.r) 
  }
  
  
  
  
  
  row.weight = apply(set.data, 1, function(x)
  {
    #  sqrt(t(x - mean(x)) %*% (x - mean(x)))
    
    sum((x-mean(x)) ^ 2)
  })
  
  
  
#   sum.sum.r = 0
#   for (i in 1:nrow(set.data)) 
#     for (j in 1:nrow(set.data))
#     {
#       x.i = set.data[i,]
#       x.j = set.data[j,]
#       
#       d.x.i = x.i - mean(x.i)
#       d.x.j = x.j - mean(x.j)
#       sum.sum.r = sum.sum.r +
#         sum(d.x.i * d.x.j) /
#         sqrt(row.weight[i] * row.weight[j])
#     }

  
  
  x.j = set.data  
  d.x.j = x.j - rowMeans(x.j)
  
  if (missing("weights"))
  {  
    sum.sum.r = sum(sapply(1:nrow(set.data), function(i)
    {
      x.i = set.data[i,]      
      d.x.i = x.i - mean(x.i)
      
      return(sum(colSums(t(d.x.j) * d.x.i) / sqrt(row.weight[i] * row.weight)))
    }))
    
  } else
  {
    
    sum.sum.r = sum(sapply(1:nrow(set.data), function(i)
    {
      x.i = set.data[i,]      
      d.x.i = x.i - mean(x.i)
      
      return(sum(weights[i]*weights*colSums(t(d.x.j) * d.x.i) / sqrt(row.weight[i] * row.weight)))
    }))
    
    sum.sum.r = 1/sum(weights)^2 * sum.sum.r
  }
  

  
  

  
  
  
  beta.score = (sum.r / sqrt(sum.sum.r))^2
  
  beta.significance = 1 - pbeta(beta.score, 0.5, (ncol(set.data)-2) / 2)
  
  
  return(list(beta.score=as.vector(beta.score), beta.significance=as.vector(beta.significance)))
  
}

  
