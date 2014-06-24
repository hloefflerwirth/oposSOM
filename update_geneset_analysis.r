
if (is.null(env$spot.list.group.overexpression)) {
  stop("This environment seems to be pretty old. It is missing \"spot.list.group.overexpression\".")
}

library(oposSOM)
options(error=quote({dump.frames(to.file=TRUE)}))

dir.create(paste(env$files.name, "- Results"), showWarnings=F)
dir.create(paste(env$files.name, "- Results/CSV Sheets"), showWarnings=F)
dir.create(paste(env$files.name, "- Results/Geneset Analysis"), showWarnings=F)
dir.create(paste(env$files.name, "- Results/3rd lvl Spot Analysis"), showWarnings=F)

env$preferences$geneset.analysis <- T
env$preferences$geneset.analysis.exact <- (exists("doAll") && doAll == TRUE)
doAll <- env$preferences$geneset.analysis.exact

if (!doAll) {
  # Backup old results
  gs.def.list <- env$gs.def.list
  spot.list.samples <- env$spot.list.samples
  spot.list.correlation <- env$spot.list.correlation
  spot.list.group.overexpression <- env$spot.list.group.overexpression
  spot.list.kmeans <- env$spot.list.kmeans
  spot.list.overexpression <- env$spot.list.overexpression
  spot.list.underexpression <- env$spot.list.underexpression
}

# Get new gene sets
oposSOM:::util.info("Preparing Annotation")
oposSOM:::util.call(oposSOM:::pipeline.prepareAnnotation, env)

if (!doAll) {
  # Extract new genesets
  env$gs.def.list <- env$gs.def.list[setdiff(names(env$gs.def.list), names(gs.def.list))]
  env$gs.def.list.categories <- sapply(gs.def.list, function(x) { x$Type })
  oposSOM:::util.info("Got", length(env$gs.def.list), "new GO sets")
}

# Calc
oposSOM:::util.call(oposSOM:::pipeline.genesetStatisticSamples, env)
oposSOM:::util.call(oposSOM:::pipeline.genesetStatisticIntegral, env)
oposSOM:::util.call(oposSOM:::pipeline.genesetProfilesAndMaps, env)

if (!doAll) {
  # Merge new results
  env$gs.def.list <- modifyList(gs.def.list, env$gs.def.list)
  env$gs.def.list.categories <- sapply(gs.def.list, function(x) { x$Type })

  for (i in seq_along(spot.list.samples)) {
    env$spot.list.samples[[i]]$GSZ.scores <-
      c(spot.list.samples[[i]]$GSZ.scores,
        env$spot.list.samples[[i]]$GSZ.scores)

    for (j in seq_along(spot.list.samples[[i]]$spots)) {
      env$spot.list.samples[[i]]$spots[[j]]$GSZ.score <-
        c(spot.list.samples[[i]]$spots[[j]]$GSZ.score,
          env$spot.list.samples[[i]]$spots[[j]]$GSZ.score)

      env$spot.list.samples[[i]]$spots[[j]]$Fisher.p <-
        c(spot.list.samples[[i]]$spots[[j]]$Fisher.p,
          env$spot.list.samples[[i]]$spots[[j]]$Fisher.p)
    }
  }

  for (i in seq_along(spot.list.correlation$spots)) {
    env$spot.list.correlation$spots[[i]]$Fisher.p <-
      c(spot.list.correlation$spots[[i]]$Fisher.p,
        env$spot.list.correlation$spots[[i]]$Fisher.p)
  }

  for (i in seq_along(spot.list.group.overexpression$spots)) {
    env$spot.list.group.overexpression$spots[[i]]$Fisher.p <-
      c(spot.list.group.overexpression$spots[[i]]$Fisher.p,
        env$spot.list.group.overexpression$spots[[i]]$Fisher.p)
  }

  for (i in seq_along(spot.list.kmeans$spots)) {
    env$spot.list.kmeans$spots[[i]]$Fisher.p <-
      c(spot.list.kmeans$spots[[i]]$Fisher.p,
        env$spot.list.kmeans$spots[[i]]$Fisher.p)
  }

  for (i in seq_along(spot.list.overexpression$spots)) {
    env$spot.list.overexpression$spots[[i]]$Fisher.p <-
      c(spot.list.overexpression$spots[[i]]$Fisher.p,
        env$spot.list.overexpression$spots[[i]]$Fisher.p)
  }

  for (i in seq_along(spot.list.underexpression$spots)) {
    env$spot.list.underexpression$spots[[i]]$Fisher.p <-
      c(spot.list.underexpression$spots[[i]]$Fisher.p,
        env$spot.list.underexpression$spots[[i]]$Fisher.p)
  }
}

# Calc GSZ scores
env$samples.GSZ.scores <- do.call(cbind, lapply(env$spot.list.samples, function(x)
{
  x$GSZ.score[names(env$gs.def.list)]
}))

# And write them to disk
filename <- file.path(paste(env$files.name, "- Results"), "CSV Sheets", "Sample GSZ scores.csv")
util.info("Writing:", filename)
write.csv2(env$samples.GSZ.scores, filename)

# Output
oposSOM:::util.call(oposSOM:::pipeline.genesetOverviews, env)
oposSOM:::util.call(oposSOM:::pipeline.geneLists, env)
oposSOM:::util.call(oposSOM:::pipeline.summarySheetsSamples, env)
oposSOM:::util.call(oposSOM:::pipeline.summarySheetsIntegral, env)
oposSOM:::util.call(oposSOM:::pipeline.3rdLvlSummarySheets, env)
oposSOM:::util.call(oposSOM:::pipeline.htmlGenesetAnalysis, env)
save(env, file=paste(env$files.name, ".RData", sep=""))

# More
oposSOM:::util.call(oposSOM:::pipeline.groupAnalysis, env)
oposSOM:::util.call(oposSOM:::pipeline.differenceAnalyses, env)
