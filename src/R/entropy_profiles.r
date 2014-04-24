pipeline.entropyProfiles <- function()
{
  filename <- file.path(paste(files.name, "- Results"), "Supporting Maps&Profiles", "Entropy Profiles.pdf")
  util.info("Writing:", filename)

  ### Metagene Mean Expresion + Variance ###
  pdf(filename, 42/2.54, 21/2.54)
  par(mar=c(10, 6, 4, 5))

  barplot(apply(metadata, 2, mean), col=group.colors,
          main="Metagene Mean Expression", names.arg=colnames(indata), las=2,
          cex.main=2.5, cex.lab=2, cex.axis=2, cex.names=1.2,
          border=if (ncol(indata) < 80) "black" else NA)

  box()

  if (length(unique(group.labels)) > 1)
  {
    mean.boxes <- by(apply(metadata, 2, mean), group.labels, c)[unique(group.labels)]

    boxplot(mean.boxes, col=unique.group.colors, las=2,
            main="Metagene Mean Expression", cex.main=2.5, cex.axis=2, xaxt="n")

    axis(1, 1:length(unique.group.colors), unique(group.labels), las=2)
  }

  barplot(apply(metadata, 2, var), col=group.colors, main="Metagene Variance",
          names.arg=colnames(indata), las=2, cex.main=2.5, cex.lab=2,
          cex.axis=2, cex.names=1.2, border=if (ncol(indata) < 80) "black" else NA)

  box()

  if (length(unique(group.labels)) > 1)
  {
    mean.sd.boxes < by(apply(metadata, 2, var), group.labels, c)[unique(group.labels)]

    boxplot(mean.sd.boxes, col=unique.group.colors, las=2,
            main="Metagene Variance", cex.main=2.5, cex.axis=2, xaxt="n")

    axis(1, 1:length(unique.group.colors), unique(group.labels), las=2)
  }

  ### Entropy ###
  q25 <- quantile(metadata,0.25)
  q75 <- quantile(metadata,0.75)

  p.metadata <- apply(metadata, 2, function(x)
  {
    hist(x, breaks=c(min(x), q25, q75, max(x)), plot=F)$counts / preferences$dim.som1^2
  })

  q25 <- quantile(metadata * som.result$code.sum[,"nobs"],0.25)
  q75 <- quantile(metadata * som.result$code.sum[,"nobs"],0.75)

  p.metadata.weighted <- apply(metadata * som.result$code.sum[,"nobs"], 2, function(x)
  {
    hist(x, breaks=c(min(x), q25, q75, max(x)), plot=F)$counts / preferences$dim.som1^2
  })

  ### Standard sample-related metagene entropy
  H <- apply(p.metadata, 2, function(p)
  {
    p <- p[p != 0]
    -sum(p * log2(p))
  })

  ylim <- c(min(H)-0.1*(max(H)-min(H)),max(H)+0.1*(max(H)-min(H)))

  barplot(H, col=group.colors, main="Standard Metagene Entropy",
          names.arg=colnames(indata), las=2, cex.main=2.5, cex.lab=2,
          cex.axis=2, cex.names=1.2, ylim=ylim, xpd=F,
          border=if (ncol(indata) < 80) "black" else NA)

  box()

  if (length(unique(group.labels)) > 1)
  {
    mean.sd.boxes <- by(H, group.labels, c)[unique(group.labels)]

    boxplot(mean.sd.boxes, col=unique.group.colors, las=2,
            main="Standard Metagene Entropy", cex.main=2.5, cex.axis=2, xaxt="n")

    axis(1, 1:length(unique.group.colors), unique(group.labels), las=2)
  }

  ### Weighted sample-related metagene entropy
  H <- apply(p.metadata.weighted, 2, function(p)
  {
    p <- p[p != 0]
    -sum(p * log2(p))
  })

  ylim <- c(min(H)-0.1*(max(H)-min(H)),max(H)+0.1*(max(H)-min(H)))

  barplot(H, col=group.colors, main="Weighted Metagene Entropy",
          names.arg=colnames(indata), las=2, cex.main=2.5, cex.lab=2,
          cex.axis=2, cex.names=1.2, ylim=ylim, xpd=F,
          border=if (ncol(indata) < 80) "black" else NA)

  box()

  if (length(unique(group.labels)) > 1)
  {
    mean.sd.boxes <- by(H, group.labels, c)[unique(group.labels)]

    boxplot(mean.sd.boxes, col=unique.group.colors, las=2,
            main="Weighted Metagene Entropy", cex.main=2.5, cex.axis=2, xaxt="n")

    axis(1, 1:length(unique.group.colors), unique(group.labels), las=2)
  }

  ### Tsallis Entropy
  q <- 0.25

  H <- apply(p.metadata, 2, function(p)
  {
    (1/(q-1)) * (1-sum(p^q))
  })

  ylim <- c(min(H)-0.1*(max(H)-min(H)),max(H)+0.1*(max(H)-min(H)))

  barplot(H, col=group.colors, main="Tsallis Metagene Entropy",
          names.arg=colnames(indata), las=2, cex.main=2.5, cex.lab=2,
          cex.axis=2, cex.names=1.2, ylim=ylim, xpd=F,
          border=if (ncol(indata) < 80) "black" else NA)

  box()

  if (length(unique(group.labels)) > 1)
  {
    mean.sd.boxes <- by(H, group.labels, c)[unique(group.labels)]

    boxplot(mean.sd.boxes, col=unique.group.colors, las=2,
            main="Tsallis Metagene Entropy", cex.main=2.5, cex.axis=2, xaxt="n")

    axis(1, 1:length(unique.group.colors), unique(group.labels), las=2)
  }

  dev.off()
}
