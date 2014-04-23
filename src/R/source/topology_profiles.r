


metadata.scaled = apply(metadata, 2, function(x){ (x-min(x))/(max(x)-min(x)) })
metadata.scaled.loglog = apply(loglog.metadata, 2, function(x){ (x-min(x))/(max(x)-min(x)) })





if (length(unique(group.labels)) > 1)
{
pdf(paste(files.name, "- Results/Supporting Maps&Profiles/Topology Profiles.pdf"), 42/2.54, 21/2.54)
  par(mar=c(10,6,4,5), mfrow=c(1,2))
} else
{
pdf(paste(files.name, "- Results/Supporting Maps&Profiles/Topology Profiles.pdf"), 21/2.54, 21/2.54)
  par(mar=c(10,6,4,5))
}





  ### Number of red spots ###


  n.spots = sapply(GS.infos.samples, function(x)
  {
    sum(sapply(x$spots, function(x){ if (x$type=="overexpressed") 1 else 0 }))
  })

  barplot(n.spots, col=group.colors, main="Number of overexpressed spots (logFC)", names.arg=colnames(indata), las=2, cex.main=2.5, cex.lab=2, cex.axis=2, cex.names=1.2, border = if (ncol(indata) < 80) "black" else NA)
    box()

  if (length(unique(group.labels)) > 1)
  {
    mean.boxes = by(n.spots, group.labels, c)[unique(group.labels)]
    boxplot(mean.boxes, col=unique.group.colors, las=2, main="Number of overexpressed spots (logFC)", cex.main=2.5, cex.axis=2, xaxt="n")
      axis(1, 1:length(unique.group.colors), unique(group.labels), las=2)
  }


  ### Spot number distribution ###

  n.spots.groups = sapply(tapply(n.spots, group.labels, c), function(x) hist(x, breaks=c(0:max(n.spots)), plot=F)$counts)
  if (is.vector(n.spots.groups))
  {
    n = names(n.spots.groups)
    n.spots.groups = matrix(n.spots.groups,nrow=1)
    colnames(n.spots.groups) = n
  }
  n.spots.groups = sapply(unique(group.labels), function(x) return(n.spots.groups[,x] / table(group.labels)[x]))

  par(mfrow=c(1,1))
  barplot(as.vector(n.spots.groups), col=rep(unique.group.colors,each=max(n.spots)), names.arg=rep(c(1:max(n.spots)),length(unique(group.labels))), las=1, main="Relative numbers of overexpressed spots (logFC)",cex.main=2.5, cex.lab=2, cex.axis=2, cex.names=0.9, border = if (ncol(indata) < 80) "black" else NA, ylim=c(0,1))
  par(mfrow=c(1,2))




  ### Fraction of red metagenes ###

  K.red = apply(metadata.scaled, 2, function(x){ length(which(x > 0.9)) }) / preferences$dim.som1^2


  barplot(K.red, col=group.colors, main="Fraction of red metagenes (logFC)", names.arg=colnames(indata), las=2, cex.main=2.5, cex.lab=2, cex.axis=2, cex.names=1.2, border = if (ncol(indata) < 80) "black" else NA)
    box()

  if (length(unique(group.labels)) > 1)
  {
    mean.boxes = by(K.red, group.labels, c)[unique(group.labels)]
    boxplot(mean.boxes, col=unique.group.colors, las=2, main="Fraction of red metagenes (logFC)", cex.main=2.5, cex.axis=2, xaxt="n")
      axis(1, 1:length(unique.group.colors), unique(group.labels), las=2)
  }






  K.red.loglog = apply(metadata.scaled.loglog, 2, function(x){ length(which(x > 0.5)) }) / preferences$dim.som1^2


  barplot(K.red.loglog, col=group.colors, main="Fraction of red metagenes (loglogFC)", names.arg=colnames(indata), las=2, cex.main=2.5, cex.lab=2, cex.axis=2, cex.names=1.2, border = if (ncol(indata) < 80) "black" else NA)
    box()

  if (length(unique(group.labels)) > 1)
  {
    mean.boxes = by(K.red.loglog, group.labels, c)[unique(group.labels)]
    boxplot(mean.boxes, col=unique.group.colors, las=2, main="Fraction of red metagenes (loglogFC)", cex.main=2.5, cex.axis=2, xaxt="n")
      axis(1, 1:length(unique.group.colors), unique(group.labels), las=2)
  }












  ### Fraction of genes in red metagenes ###


  f = apply(metadata.scaled, 2, function(x){ sum(som.result$code.sum[which(x > 0.9), "nobs"]) }) / nrow(indata)



  barplot(f, col=group.colors, main="Fraction of red genes (logFC)", names.arg=colnames(indata), las=2, cex.main=2.5, cex.lab=2, cex.axis=2, cex.names=1.2, border = if (ncol(indata) < 80) "black" else NA)
    box()

  if (length(unique(group.labels)) > 1)
  {
    mean.boxes = by(f, group.labels, c)[unique(group.labels)]
    boxplot(mean.boxes, col=unique.group.colors, las=2, main="Fraction of red genes (logFC)", cex.main=2.5, cex.axis=2, xaxt="n")
      axis(1, 1:length(unique.group.colors), unique(group.labels), las=2)
  }







  f = apply(metadata.scaled.loglog, 2, function(x){ sum(som.result$code.sum[which(x > 0.5), "nobs"]) }) / nrow(indata)



  barplot(f, col=group.colors, main="Fraction of red genes (loglogFC)", names.arg=colnames(indata), las=2, cex.main=2.5, cex.lab=2, cex.axis=2, cex.names=1.2, border = if (ncol(indata) < 80) "black" else NA)
    box()

  if (length(unique(group.labels)) > 1)
  {
    mean.boxes = by(f, group.labels, c)[unique(group.labels)]
    boxplot(mean.boxes, col=unique.group.colors, las=2, main="Fraction of red genes (loglogFC)", cex.main=2.5, cex.axis=2, xaxt="n")
      axis(1, 1:length(unique.group.colors), unique(group.labels), las=2)
  }








  ### Length of borderline ###



  K.border = apply(metadata.scaled, 2, function(x)
  {
    m = matrix(x, preferences$dim.som1, preferences$dim.som1)
    m[which(m < 0.9)] = NA

    count.border.metagenes = 0
    for (i in 1:preferences$dim.som1)
      for (j in 1:preferences$dim.som1)
        if (!is.na(m[i,j]))
        {
          neighbours = sapply(get.neighbors(i, j, preferences$dim.som1), function(x){ m[x[1],x[2]] })

          if (length(which(is.na(neighbours))))
          {
            count.border.metagenes = count.border.metagenes + 1
          } else
          if (i %in% c(1,preferences$dim.som1) || j %in% c(1,preferences$dim.som1))
          {
            count.border.metagenes = count.border.metagenes + 1
          }
        }


    return(count.border.metagenes)

  })


  barplot(K.border, col=group.colors, main="Length of borderline (logFC)", names.arg=colnames(indata), las=2, cex.main=2.5, cex.lab=2, cex.axis=2, cex.names=1.2, border = if (ncol(indata) < 80) "black" else NA)
    box()

  if (length(unique(group.labels)) > 1)
  {
    mean.boxes = by(K.border, group.labels, c)[unique(group.labels)]
    boxplot(mean.boxes, col=unique.group.colors, las=2, main="Length of borderline (logFC)", cex.main=2.5, cex.axis=2, xaxt="n")
      axis(1, 1:length(unique.group.colors), unique(group.labels), las=2)
  }









  K.border.loglog = apply(metadata.scaled.loglog, 2, function(x)
  {
    m = matrix(x, preferences$dim.som1, preferences$dim.som1)
    m[which(m < 0.9)] = NA

    count.border.metagenes = 0
    for (i in 1:preferences$dim.som1)
      for (j in 1:preferences$dim.som1)
        if (!is.na(m[i,j]))
        {
          neighbours = sapply(get.neighbors(i, j, preferences$dim.som1), function(x){ m[x[1],x[2]] })

          if (length(which(is.na(neighbours))))
          {
            count.border.metagenes = count.border.metagenes + 1
          } else
          if (i %in% c(1,preferences$dim.som1) || j %in% c(1,preferences$dim.som1))
          {
            count.border.metagenes = count.border.metagenes + 1
          }
        }


    return(count.border.metagenes)

  })


  barplot(K.border.loglog, col=group.colors, main="Length of borderline (loglogFC)", names.arg=colnames(indata), las=2, cex.main=2.5, cex.lab=2, cex.axis=2, cex.names=1.2, border = if (ncol(indata) < 80) "black" else NA)
    box()

  if (length(unique(group.labels)) > 1)
  {
    mean.boxes = by(K.border.loglog, group.labels, c)[unique(group.labels)]
    boxplot(mean.boxes, col=unique.group.colors, las=2, main="Length of borderline (loglogFC)", cex.main=2.5, cex.axis=2, xaxt="n")
      axis(1, 1:length(unique.group.colors), unique(group.labels), las=2)
  }















  ### Compactness of spots ###


  C = K.red / K.border


  barplot(C, col=group.colors, main="Compactness of spots (logFC)", names.arg=colnames(indata), las=2, cex.main=2.5, cex.lab=2, cex.axis=2, cex.names=1.2, border = if (ncol(indata) < 80) "black" else NA)
    box()

  if (length(unique(group.labels)) > 1)
  {
    mean.boxes = by(C, group.labels, c)[unique(group.labels)]
    boxplot(mean.boxes, col=unique.group.colors, las=2, main="Compactness of spots (logFC)", cex.main=2.5, cex.axis=2, xaxt="n")
      axis(1, 1:length(unique.group.colors), unique(group.labels), las=2)
  }




  C = K.red.loglog / K.border.loglog


  barplot(C, col=group.colors, main="Compactness of spots (loglogFC)", names.arg=colnames(indata), las=2, cex.main=2.5, cex.lab=2, cex.axis=2, cex.names=1.2, border = if (ncol(indata) < 80) "black" else NA)
    box()

  if (length(unique(group.labels)) > 1)
  {
    mean.boxes = by(C, group.labels, c)[unique(group.labels)]
    boxplot(mean.boxes, col=unique.group.colors, las=2, main="Compactness of spots (loglogFC)", cex.main=2.5, cex.axis=2, xaxt="n")
      axis(1, 1:length(unique.group.colors), unique(group.labels), las=2)
  }








  ### Shape of spots ###


  C = (K.red * preferences$dim.som1^2) / K.border^2


  barplot(C, col=group.colors, main="Shape of spots (logFC)", names.arg=colnames(indata), las=2, cex.main=2.5, cex.lab=2, cex.axis=2, cex.names=1.2, border = if (ncol(indata) < 80) "black" else NA)
    box()

  if (length(unique(group.labels)) > 1)
  {
    mean.boxes = by(C, group.labels, c)[unique(group.labels)]
    boxplot(mean.boxes, col=unique.group.colors, las=2, main="Shape of spots (logFC)", cex.main=2.5, cex.axis=2, xaxt="n")
      axis(1, 1:length(unique.group.colors), unique(group.labels), las=2)
  }




  C = (K.red.loglog * preferences$dim.som1^2) / K.border.loglog^2


  barplot(C, col=group.colors, main="Shape of spots (loglogFC)", names.arg=colnames(indata), las=2, cex.main=2.5, cex.lab=2, cex.axis=2, cex.names=1.2, border = if (ncol(indata) < 80) "black" else NA)
    box()

  if (length(unique(group.labels)) > 1)
  {
    mean.boxes = by(C, group.labels, c)[unique(group.labels)]
    boxplot(mean.boxes, col=unique.group.colors, las=2, main="Shape of spots (loglogFC)", cex.main=2.5, cex.axis=2, xaxt="n")
      axis(1, 1:length(unique.group.colors), unique(group.labels), las=2)
  }




dev.off()






































