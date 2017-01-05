pipeline.topologyProfiles <- function()
{
  metadata.scaled <- apply(metadata, 2, function(x) { (x-min(x)) / (max(x)-min(x)) })

  filename <- file.path(paste(files.name, "- Results"), "Supporting Maps&Profiles", "Topology Profiles.pdf")
  util.info("Writing:", filename)

  if (length(unique(group.labels)) > 1)
  {
    pdf(filename, 42/2.54, 21/2.54)
    par(mar=c(10, 6, 4, 5), mfrow=c(1, 2))
  } else
  {
    pdf(filename, 21/2.54, 21/2.54)
    par(mar=c(10, 6, 4, 5))
  }

  ### Number of overexpression spots ###
  spotdata <- get(paste("spot.list.",preferences$standard.spot.modules,sep=""))$spotdata
  n.spots <- colSums( spotdata>sd(spotdata) )

  barplot(n.spots, col=group.colors, main="Number of activated spot modules",
          names.arg="", las=2, cex.main=2, 
          border=if (ncol(indata) < 80) "black" else NA)
  box()

  if (length(unique(group.labels)) > 1)
  {
    mean.boxes <- by(n.spots, group.labels, c)[unique(group.labels)]

    boxplot(mean.boxes, col=groupwise.group.colors, las=2, xaxt="n")

    axis(1, seq_along(groupwise.group.colors), unique(group.labels), las=2)
  }

  ### Spot number distribution ###
  n.spots.groups <- sapply(tapply(n.spots, group.labels, c), function(x)
  {
    ret <- table(x)[as.character(0:max(n.spots))]
    ret[which(is.na(ret))] <- 0
    names(ret) <- as.character(0:max(n.spots))
    return(ret)
  })

  if (is.vector(n.spots.groups))
  {
    n.spots.groups <- matrix(n.spots.groups,nrow=1)
    colnames(n.spots.groups) <- unique(group.labels)
  }

  n.spots.groups <- sapply(unique(group.labels), function(x)
  {
    return(n.spots.groups[,x] / table(group.labels)[x])
  })

  par(mfrow=c(1, 1))

  barplot(as.vector(n.spots.groups), col=rep(groupwise.group.colors,each=max(n.spots)+1),
          names.arg=rep(c(0:max(n.spots)), length(unique(group.labels))), las=1,
          main="Fraction of samples showing respective number of overexpression spots",
          cex.main=2, border=if (ncol(indata) < 80) "black" else NA, ylim=c(0, 1))

  par(mfrow=c(1, 2))

  ### Fraction of red metagenes ###
  K.red <- apply(metadata.scaled, 2, function(x) { length(which(x > 0.9)) }) / preferences$dim.1stLvlSom^2

  barplot(K.red, col=group.colors, main="Fraction of overexpressed metagenes", cex.main=2,
          names.arg="", las=2, border=if (ncol(indata) < 80) "black" else NA)

  box()

  if (length(unique(group.labels)) > 1)
  {
    mean.boxes <- by(K.red, group.labels, c)[unique(group.labels)]

    boxplot(mean.boxes, col=groupwise.group.colors, las=2, main="", xaxt="n")

    axis(1, seq_along(groupwise.group.colors), unique(group.labels), las=2)
  }

  ### Length of borderline ###
  K.border <- apply(metadata.scaled, 2, function(x)
  {
    m <- matrix(x, preferences$dim.1stLvlSom, preferences$dim.1stLvlSom)
    m[which(m < 0.9)] <- NA
    count.border.metagenes <- 0

    for (i in 1:preferences$dim.1stLvlSom)
    {
      for (j in 1:preferences$dim.1stLvlSom)
      {
        if (!is.na(m[i,j]))
        {
          neighbours <- sapply(get.neighbors(i, j, preferences$dim.1stLvlSom),
                               function(x) { m[x[1],x[2]] })

          if (length(which(is.na(neighbours))))
          {
            count.border.metagenes <- count.border.metagenes + 1
          } else if (i %in% c(1, preferences$dim.1stLvlSom) || j %in% c(1, preferences$dim.1stLvlSom))
          {
            count.border.metagenes <- count.border.metagenes + 1
          }
        }
      }
    }

    return(count.border.metagenes)
  })

  barplot(K.border, col=group.colors, main="Length of borderline along overexpressed metagenes",
          names.arg="", las=2, cex.main=1.8, 
          border=if (ncol(indata) < 80) "black" else NA)

  box()

  if (length(unique(group.labels)) > 1)
  {
    mean.boxes = by(K.border, group.labels, c)[unique(group.labels)]

    boxplot(mean.boxes, col=groupwise.group.colors, las=2, main="", xaxt="n")

    axis(1, seq_along(groupwise.group.colors), unique(group.labels), las=2)
  }

  ### Compactness of spots ###
  C <- K.red / K.border

  barplot(C, col=group.colors, main="Compactness of spots",
          names.arg="", las=2, cex.main=2, 
          border=if (ncol(indata) < 80) "black" else NA)

  box()

  if (length(unique(group.labels)) > 1)
  {
    mean.boxes <- by(C, group.labels, c)[unique(group.labels)]

    boxplot(mean.boxes, col=groupwise.group.colors, las=2, main="", xaxt="n")

    axis(1, seq_along(groupwise.group.colors), unique(group.labels), las=2)
  }

  ### Shape of spots ###
  C <- (K.red * preferences$dim.1stLvlSom^2) / K.border^2

  barplot(C, col=group.colors, main="Shape of spots",
          names.arg="", las=2, cex.main=2, 
          border=if (ncol(indata) < 80) "black" else NA)

  box()

  if (length(unique(group.labels)) > 1)
  {
    mean.boxes <- by(C, group.labels, c)[unique(group.labels)]

    boxplot(mean.boxes, col=groupwise.group.colors, las=2, main="",xaxt="n")

    axis(1, seq_along(groupwise.group.colors), unique(group.labels), las=2)
  }

  dev.off()
}
