  

if (length(unique(group.labels)) > 1)
{
    
    group.assignment = match(group.labels, unique(group.labels))
    names(group.assignment) = colnames(indata)
  
    
     bootstrap.error = rep(0,ncol(indata))
     names(bootstrap.error) = colnames(indata)
    
    n.bootstrap = min(1000, 5*ncol(indata))
     for (i in 1:n.bootstrap)
     {
        resample = sample(colnames(indata), ncol(indata), replace=T)
       
        km = c()
         suppressWarnings({ km = kmeans(t(metadata[,resample]), centers=t(group.metadata), iter.max=10, algorithm = "Lloyd")  })
       
        if (length(km) > 0)
         {
           faults = which(km$cluster[resample] != group.assignment[resample])
       
          if (length(faults) > 0) bootstrap.error[names(faults)] = bootstrap.error[names(faults)] + 1     
        }
     }
    group.bootstrap.score = round(100*(1 - bootstrap.error / n.bootstrap), 1)

} else
{
  
  group.bootstrap.score = rep(100, ncol(indata))
  
}
    
  
  
    
    
  
    
    
    
    
    
    
    
    
  
  
  
  
  
  
  
  
  
  
  
  
