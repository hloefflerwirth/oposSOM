pipeline.entropyProfiles <- function(env)
{
  filename <- file.path("Supporting Maps&Profiles", "Entropy Profiles.pdf")
  util.info("Writing:", filename)

  ### Metagene Mean Expresion + Variance ###
  pdf(filename, 42/2.54, 21/2.54, useDingbats=FALSE)
  par(mar=c(10, 6, 4, 5))

  barplot(apply(env$metadata, 2, mean), col=env$group.colors,
          main=bquote("Mean metagene expression: <e"[m]^meta~">"), names.arg="", las=2,
          cex.main=2.5,border=if (ncol(env$indata) < 80) "black" else NA)
  box()

  if (length(unique(env$group.labels)) > 1)
  {
    mean.boxes <- by(apply(env$metadata, 2, mean), env$group.labels, c)[unique(env$group.labels)]

    boxplot(mean.boxes, col=env$groupwise.group.colors, las=2,
            main=bquote("Mean metagene expression: <e"[m]^meta~">"), cex.main=2.5, xaxt="n")
  
    axis(1, seq_along(env$groupwise.group.colors), unique(env$group.labels), las=2)
  }

  barplot(apply(env$metadata, 2, var), col=env$group.colors, main=bquote("Metagene variance: var(e"[m]^meta~")"),
          names.arg="", las=2, cex.main=2.5, border=if (ncol(env$indata) < 80) "black" else NA)
  box()

  if (length(unique(env$group.labels)) > 1)
  {
    mean.sd.boxes <- by(apply(env$metadata, 2, var), env$group.labels, c)[unique(env$group.labels)]

    boxplot(mean.sd.boxes, col=env$groupwise.group.colors, las=2,
            main=bquote("Metagene variance: var(e"[m]^meta~")"), cex.main=2.5, xaxt="n")

    axis(1, seq_along(env$groupwise.group.colors), unique(env$group.labels), las=2)
  }

  ### Entropy ###
  q25 <- quantile(env$metadata,0.25)
  q75 <- quantile(env$metadata,0.75)

  p.metadata <- apply(env$metadata, 2, function(x)
  {
    hist(x, breaks=c(min(x), q25, q75, max(x)), plot=FALSE)$counts / env$preferences$dim.1stLvlSom^2
  })

  q25 <- quantile(env$metadata * env$som.result$node.summary[,"n.features"],0.25)
  q75 <- quantile(env$metadata * env$som.result$node.summary[,"n.features"],0.75)

  p.metadata.weighted <- apply(env$metadata * env$som.result$node.summary[,"n.features"], 2, function(x)
  {
    hist(x, breaks=c(min(x), q25, q75, max(x)), plot=FALSE)$counts / env$preferences$dim.1stLvlSom^2
  })

  ### Standard sample-related metagene entropy
  H <- apply(p.metadata, 2, function(p)
  {
    p <- p[p != 0]
    -sum(p * log2(p))
  })

  ylim <- c(min(H)-0.1*(max(H)-min(H)),max(H)+0.1*(max(H)-min(H)))

  barplot(H, col=env$group.colors, main=bquote("Standard metagene entropy: h"[m]),
          names.arg="", las=2, cex.main=2.5,
          ylim=ylim, xpd=FALSE, border=if (ncol(env$indata) < 80) "black" else NA)

  box()

  if (length(unique(env$group.labels)) > 1)
  {
    mean.sd.boxes <- by(H, env$group.labels, c)[unique(env$group.labels)]

    boxplot(mean.sd.boxes, col=env$groupwise.group.colors, las=2,
            main=bquote("Standard metagene entropy: h"[m]~""), cex.main=2.5, xaxt="n")

    axis(1, seq_along(env$groupwise.group.colors), unique(env$group.labels), las=2)
  }

  ### Weighted sample-related metagene entropy
  H <- apply(p.metadata.weighted, 2, function(p)
  {
    p <- p[p != 0]
    -sum(p * log2(p))
  })

  ylim <- c(min(H)-0.1*(max(H)-min(H)),max(H)+0.1*(max(H)-min(H)))

  barplot(H, col=env$group.colors, main=bquote("Metagene entropy weighted for metagene population: h"[m]^weighted),
          names.arg="", las=2, cex.main=2.5, ylim=ylim, 
          xpd=FALSE, border=if (ncol(env$indata) < 80) "black" else NA)

  box()

  if (length(unique(env$group.labels)) > 1)
  {
    mean.sd.boxes <- by(H, env$group.labels, c)[unique(env$group.labels)]

    boxplot(mean.sd.boxes, col=env$groupwise.group.colors, las=2,
            main=bquote("Metagene entropy weighted for metagene population: h"[m]^weighted~""), cex.main=2.5, xaxt="n")

    axis(1, seq_along(env$groupwise.group.colors), unique(env$group.labels), las=2)
  }

  ### Tsallis Entropy
  q <- 0.25

  H <- apply(p.metadata, 2, function(p)
  {
    (1/(q-1)) * (1-sum(p^q))
  })

  ylim <- c(min(H)-0.1*(max(H)-min(H)),max(H)+0.1*(max(H)-min(H)))

  barplot(H, col=env$group.colors, main=bquote("Tsallis metagene entropy: h"[m]^Tsallis~""),
          names.arg="", las=2, cex.main=2.5, ylim=ylim, 
          xpd=FALSE, border=if (ncol(env$indata) < 80) "black" else NA)

  box()

  if (length(unique(env$group.labels)) > 1)
  {
    mean.sd.boxes <- by(H, env$group.labels, c)[unique(env$group.labels)]

    boxplot(mean.sd.boxes, col=env$groupwise.group.colors, las=2,
            main=bquote("Tsallis metagene entropy: h"[m]^Tsallis~""), cex.main=2.5, xaxt="n")

    axis(1, seq_along(env$groupwise.group.colors), unique(env$group.labels), las=2)
  }

  dev.off()
}
