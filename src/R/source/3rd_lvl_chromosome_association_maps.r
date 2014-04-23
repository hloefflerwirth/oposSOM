

  
  all.chromos = sort.label(names(gene.positions.list))
  all.chromos = all.chromos[which(nchar(all.chromos) < 4)]
   
  all.chromos.pq = paste(unlist(sapply(all.chromos, function(x) rep(x, length(gene.positions.list[[x]])))), unlist(sapply(gene.positions.list[all.chromos],names)))
  all.chromos.pq = unique(sub("[0-9]+$", "", all.chromos.pq))
  
  all.chromos.pq.count = table(sub("[0-9]+$", "", na.omit(gene.positions)))[all.chromos.pq]
  
  spot.gene.chromo.pq.matrix = t(sapply(set.list$spots, function(x) 
  {
    spot.genes = x$genes
    spot.gene.chromo = sub("[0-9]+$", "", na.omit(gene.positions[spot.genes]))
    spot.gene.chromo = spot.gene.chromo[which(nchar(spot.gene.chromo) <= 4)]
    
    spot.gene.chromo.tab = table(spot.gene.chromo)[all.chromos.pq]
    names(spot.gene.chromo.tab) = all.chromos.pq
    spot.gene.chromo.tab[which(is.na(spot.gene.chromo.tab))] = 0
      
    return(spot.gene.chromo.tab)
  }))
  
  
  spot.gene.chromo.pq.matrix = spot.gene.chromo.pq.matrix[, which(colSums(spot.gene.chromo.pq.matrix) > 0)]
  
  spot.gene.chromo.pq.matrix.scaled = t(apply(spot.gene.chromo.pq.matrix, 1, function(x) x/all.chromos.pq.count[names(x)]))
  
  
  
  
  
  
  
  max.enrichment = max(spot.gene.chromo.pq.matrix^2)
  max.enrichment.scaled = max(spot.gene.chromo.pq.matrix.scaled^2)
  
  
  spot.gene.chromo.pq.matrix = rbind(spot.gene.chromo.pq.matrix, colSums(spot.gene.chromo.pq.matrix) / max(colSums(spot.gene.chromo.pq.matrix)) * sqrt(max.enrichment))
  rownames(spot.gene.chromo.pq.matrix)[nrow(spot.gene.chromo.pq.matrix)] = "summary"
  
  spot.gene.chromo.pq.matrix.scaled = rbind(spot.gene.chromo.pq.matrix.scaled, colSums(spot.gene.chromo.pq.matrix.scaled) / max(colSums(spot.gene.chromo.pq.matrix.scaled)) * sqrt(max.enrichment.scaled))
  rownames(spot.gene.chromo.pq.matrix.scaled)[nrow(spot.gene.chromo.pq.matrix.scaled)] = "summary"
  
  
  
  # report output  

  spot.i=1
  for (spot.i in 1:nrow(spot.gene.chromo.pq.matrix))
  {

    layout(matrix(c(1,3,4,6,2,3,5,7),4), heights=c(4,1.7,12,12))
      
    
    # spot summary
    
    par(mar=c(0,0,0,0))
    plot(0, type="n", axes=F, xlab="", ylab="", xlim=c(0,1), ylim=c(0,1))

      text(0.01, 0.92, paste("Spot", rownames(spot.gene.chromo.pq.matrix)[spot.i]) , cex=2.6, adj=0)  
    
    
    if (spot.i <= length(set.list$spots))
    {
          
      samples.with.spot = names(which(set.list$spotdata[spot.i,] > sd(as.vector(set.list$spotdata))))
      samples.with.spot = table(group.labels[samples.with.spot])[unique(group.labels)]
      samples.with.spot = samples.with.spot[which(!is.na(samples.with.spot))]
      
      text(0.01, 0.7, paste("# genes in spot =", length(set.list$spots[[spot.i]]$genes), " (", sum(spot.gene.chromo.pq.matrix[spot.i,]),"with loci information)"), adj=0)
      text(0.01, 0.58, paste("# samples with spot =", sum(samples.with.spot), "(", round(100 * sum(samples.with.spot)/ncol(indata), 1), "%)"), adj=0)
      
      if (length(unique(group.labels)) > 1 && length(samples.with.spot) > 0)
      {
        for (g in 1:length(samples.with.spot))
        {                  
          text(0.02, 0.58-g*0.1, paste(names(samples.with.spot)[g], ":", samples.with.spot[g], "(", round(100 * samples.with.spot[g]/sum(group.labels == names(samples.with.spot)[g]), 1), "%)"), adj=0,  col = group.colors[match(names(samples.with.spot)[g], group.labels)])
        }
      }
    
    }
  
    
    # spot minimap
    
    if (spot.i <= length(set.list$spots))
    {
    
      par(mar=c(2,4,2,18))
      
      image(matrix(set.list$spots[[spot.i]]$mask, preferences$dim.som1, preferences$dim.som1), axes=F, col ="darkgreen")
        box()
      
    }      
    else frame()
  

    
    
    
    
    # chromosome occupation stripes
  
    if (spot.i <= length(set.list$spots))
    {
        
      par(mar=c(1,2.1,1,2.1))
      
      image(x=1:nrow(intersect.counts), y=1, z=matrix(intersect.counts[,spot.i],ncol=1), zlim=c(0,1), col=colkey, axes=F, xlab="",ylab="")
        box()
          
        abline(v=nrow(intersect.counts)-chr.sep.level+0.5, lwd=0.1, lty=3)
        axis(1, nrow(intersect.counts)-chr.lab.level,  rev(unique(chr)), las=2, tick=F, line=-0.5)
      
    }
    else frame()
    
    


    
    # chromosome association map - absolute
      
    chromo.data = spot.gene.chromo.pq.matrix
  
    par(mar=c(2.2,2.2,4,1))
    plot(0, type="n", xlab="", ylab="", xlim=c(1,ncol(chromo.data)), ylim=c(1,ncol(chromo.data)), axes=F)

    
    
      title(main="Absolute numbers", line=0.5, cex.main=1)
      box()
      
      axis(1, c(1:ncol(chromo.data)), colnames(chromo.data), las=2, cex.axis=0.7)
      axis(2, c(1:ncol(chromo.data)), colnames(chromo.data), las=2, cex.axis=0.7)
      
    
    chr.i=1
    chr.j=2
    for (chr.i in 1:ncol(chromo.data))
      for (chr.j in 1:ncol(chromo.data))
      {
        n1 = chromo.data[spot.i,chr.i]
        n2 = chromo.data[spot.i,chr.j]
        
        points(chr.i, chr.j, pch=16, cex= 0+3* n1*n2/max(chromo.data[spot.i,]^2), col="blue")
        
      }
    
  
    
    # ranked chromo list
    
    if (spot.i <= length(set.list$spots))
    {
      o = names(sort(chromo.data[spot.i,], decreasing=T))
            
      p.values=sapply(o, function(loci.i)
      {
        n.chr.spot = spot.gene.chromo.pq.matrix[spot.i,loci.i]
        n.chr.Nspot =  all.chromos.pq.count[loci.i] - n.chr.spot
        n.Nchr.spot = length(set.list$spots[[spot.i]]$genes) - n.chr.spot
        n.Nchr.Nspot = sum(all.chromos.pq.count) - n.chr.spot - n.chr.Nspot - n.Nchr.spot
        
        p = fisher.test(matrix(c(n.chr.spot,n.chr.Nspot,n.Nchr.spot,n.Nchr.Nspot),2), alternative="greater")$p.value  
        p = min(1, p)
        p = max(1e-99, p)
        
        return(p)    
      })
            
      
      x.coords = c(0.1, 0.2, 0.35, 0.5, 0.65)
      y.coords = seq(0.91, 0.04, length.out=ncol(chromo.data))
  
      par(mar=c(0,0,0,0))
      plot(0, type="n", axes=F, xlab="", ylab="", xlim=c(0,1), ylim=c(0,1))
          
        text(x.coords, rep(0.95, 5), c("Rank", "locus", "# genes", "%", "p-value"), cex=1, adj=0)
       
         text(x.coords[1], y.coords, c(1:length(o)), adj=0)    
         text(x.coords[2], y.coords, o, adj=0)           
        text(x.coords[3], y.coords, spot.gene.chromo.pq.matrix[spot.i,o], adj=0)
        text(x.coords[4], y.coords, round(100*spot.gene.chromo.pq.matrix.scaled[spot.i,o],1), adj=0)
        text(x.coords[5], y.coords, sapply(p.values,format.pval,1), adj=0)
      
      
    }
    else frame()
     
    
    
    
    
      
    
    
    
    
    # chromosome association map - relative
    
    chromo.data = spot.gene.chromo.pq.matrix.scaled  
    
    par(mar=c(2.2,2.2,4,1))
    plot(0, type="n", xlab="", ylab="", xlim=c(1,ncol(chromo.data)), ylim=c(1,ncol(chromo.data)), axes=F)  
      title(main="Relative numbers", line=0.5, cex.main=1)
      box()
      
      axis(1, c(1:ncol(chromo.data)), colnames(chromo.data), las=2, cex.axis=0.7)
      axis(2, c(1:ncol(chromo.data)), colnames(chromo.data), las=2, cex.axis=0.7)
      
    
    chr.i=1
    chr.j=2
    for (chr.i in 1:ncol(chromo.data))
      for (chr.j in 1:ncol(chromo.data))
      {
        n1 = chromo.data[spot.i,chr.i]
        n2 = chromo.data[spot.i,chr.j]
        
        points(chr.i, chr.j, pch=16, cex= 0+3* n1*n2/max(chromo.data[spot.i,]^2), col="blue")
        
      }    
    
    
    
    # ranked chromo list
    
    if (spot.i <= length(set.list$spots))
    {
      o = names(sort(chromo.data[spot.i,], decreasing=T))
      
      p.values=sapply(o, function(loci.i)
      {
        n.chr.spot = spot.gene.chromo.pq.matrix[spot.i,loci.i]
        n.chr.Nspot =  all.chromos.pq.count[loci.i] - n.chr.spot
        n.Nchr.spot = length(set.list$spots[[spot.i]]$genes) - n.chr.spot
        n.Nchr.Nspot = sum(all.chromos.pq.count) - n.chr.spot - n.chr.Nspot - n.Nchr.spot
        
        p = fisher.test(matrix(c(n.chr.spot,n.chr.Nspot,n.Nchr.spot,n.Nchr.Nspot),2), alternative="greater")$p.value  
        p = min(1, p)
        p = max(1e-99, p)
        
        return(p)    
      })
      
      
      x.coords = c(0.1, 0.2, 0.35, 0.5, 0.65)
      y.coords = seq(0.91, 0.04, length.out=ncol(chromo.data))
      
      par(mar=c(0,0,0,0))
      plot(0, type="n", axes=F, xlab="", ylab="", xlim=c(0,1), ylim=c(0,1))
      
      text(x.coords, rep(0.95, 5), c("Rank", "locus", "# genes", "%", "p-value"), cex=1, adj=0)
      
      text(x.coords[1], y.coords, c(1:length(o)), adj=0)    
      text(x.coords[2], y.coords, o, adj=0)           
      text(x.coords[3], y.coords, spot.gene.chromo.pq.matrix[spot.i,o], adj=0)
      text(x.coords[4], y.coords, round(100*spot.gene.chromo.pq.matrix.scaled[spot.i,o],1), adj=0)
      text(x.coords[5], y.coords, sapply(p.values,format.pval,1), adj=0)
      
      
    }
    else frame()
    
    
    
    
    
  }  
    

  
  





