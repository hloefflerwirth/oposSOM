pipeline.topologyProfiles <- function(env)
{
  metadata.scaled <- apply(env$metadata, 2, function(x) { (x-min(x)) / (max(x)-min(x)) })

  filename <- file.path("Supporting Maps&Profiles", "Topology Profiles.pdf")
  util.info("Writing:", filename)

  if (length(unique(env$group.labels)) > 1)
  {
    pdf(filename, 42/2.54, 21/2.54, useDingbats=FALSE)
    par(mar=c(10, 6, 4, 5), mfrow=c(1, 2))
  } else
  {
    pdf(filename, 21/2.54, 21/2.54, useDingbats=FALSE)
    par(mar=c(10, 6, 4, 5))
  }

  ### Number of overexpression spots ###
  spotdata <- env[[paste("spot.list.",env$preferences$standard.spot.modules,sep="")]]$spotdata
  n.spots <- colSums( spotdata>sd(spotdata) )

  barplot(n.spots, col=env$group.colors, main="Number of activated spot modules",
          names.arg="", las=2, cex.main=2, 
          border=if (ncol(env$indata) < 80) "black" else NA)
  box()

  if (length(unique(env$group.labels)) > 1)
  {
    mean.boxes <- by(n.spots, env$group.labels, c)[unique(env$group.labels)]

    boxplot(mean.boxes, col=env$groupwise.group.colors, las=2, xaxt="n")

    axis(1, seq_along(env$groupwise.group.colors), unique(env$group.labels), las=2)
  }

  ### Spot number distribution ###
  n.spots.groups <- sapply(tapply(n.spots, env$group.labels, c), function(x)
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

  n.spots.groups <- sapply(unique(env$group.labels), function(x)
  {
    return(n.spots.groups[,x] / table(env$group.labels)[x])
  })

  par(mfrow=c(1, 1))

  barplot(as.vector(n.spots.groups), col=rep(env$groupwise.group.colors,each=max(n.spots)+1),
          names.arg=rep(c(0:max(n.spots)), length(unique(env$group.labels))), las=1,
          main="Fraction of samples showing respective number of overexpression spots",
          cex.main=2, border=if (ncol(env$indata) < 80) "black" else NA, ylim=c(0, 1))

  par(mfrow=c(1, 2))

  ### Fraction of red metagenes ###
  K.red <- apply(metadata.scaled, 2, function(x) { length(which(x > 0.9)) }) / env$preferences$dim.1stLvlSom^2

  barplot(K.red, col=env$group.colors, main="Fraction of overexpressed metagenes", cex.main=2,
          names.arg="", las=2, border=if (ncol(env$indata) < 80) "black" else NA)

  box()

  if (length(unique(env$group.labels)) > 1)
  {
    mean.boxes <- by(K.red, env$group.labels, c)[unique(env$group.labels)]

    boxplot(mean.boxes, col=env$groupwise.group.colors, las=2, main="", xaxt="n")

    axis(1, seq_along(env$groupwise.group.colors), unique(env$group.labels), las=2)
  }

  ### Length of borderline ###
  K.border <- apply(metadata.scaled, 2, function(x)
  {
    m <- matrix(x, env$preferences$dim.1stLvlSom, env$preferences$dim.1stLvlSom)
    m[which(m < 0.9)] <- NA
    count.border.metagenes <- 0

    for (i in 1:env$preferences$dim.1stLvlSom)
    {
      for (j in 1:env$preferences$dim.1stLvlSom)
      {
        if (!is.na(m[i,j]))
        {
          neighbours <- sapply(get.neighbors(i, j, env$preferences$dim.1stLvlSom),
                               function(x) { m[x[1],x[2]] })

          if (length(which(is.na(neighbours))))
          {
            count.border.metagenes <- count.border.metagenes + 1
          } else if (i %in% c(1, env$preferences$dim.1stLvlSom) || j %in% c(1, env$preferences$dim.1stLvlSom))
          {
            count.border.metagenes <- count.border.metagenes + 1
          }
        }
      }
    }

    return(count.border.metagenes)
  })

  barplot(K.border, col=env$group.colors, main="Length of borderline along overexpressed metagenes",
          names.arg="", las=2, cex.main=1.8, 
          border=if (ncol(env$indata) < 80) "black" else NA)

  box()

  if (length(unique(env$group.labels)) > 1)
  {
    mean.boxes = by(K.border, env$group.labels, c)[unique(env$group.labels)]

    boxplot(mean.boxes, col=env$groupwise.group.colors, las=2, main="", xaxt="n")

    axis(1, seq_along(env$groupwise.group.colors), unique(env$group.labels), las=2)
  }

  ### Compactness of spots ###
  C <- K.red / K.border

  barplot(C, col=env$group.colors, main="Compactness of spots",
          names.arg="", las=2, cex.main=2, 
          border=if (ncol(env$indata) < 80) "black" else NA)

  box()

  if (length(unique(env$group.labels)) > 1)
  {
    mean.boxes <- by(C, env$group.labels, c)[unique(env$group.labels)]

    boxplot(mean.boxes, col=env$groupwise.group.colors, las=2, main="", xaxt="n")

    axis(1, seq_along(env$groupwise.group.colors), unique(env$group.labels), las=2)
  }

  ### Shape of spots ###
  C <- (K.red * env$preferences$dim.1stLvlSom^2) / K.border^2

  barplot(C, col=env$group.colors, main="Shape of spots",
          names.arg="", las=2, cex.main=2, 
          border=if (ncol(env$indata) < 80) "black" else NA)

  box()

  if (length(unique(env$group.labels)) > 1)
  {
    mean.boxes <- by(C, env$group.labels, c)[unique(env$group.labels)]

    boxplot(mean.boxes, col=env$groupwise.group.colors, las=2, main="",xaxt="n")

    axis(1, seq_along(env$groupwise.group.colors), unique(env$group.labels), las=2)
  }

  dev.off()
}
