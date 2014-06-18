library(oposSOM)

migrate <- function(fin, fout) {
  old.env <- new.env()
  load(fin, envir=old.env)

  globals <- list(GS.infos.correlation="spot.list.correlation",
                  GS.infos.group.overexpression="spot.list.group.overexpression",
                  GS.infos.kmeans="spot.list.kmeans",
                  GS.infos.overexpression="spot.list.overexpression",
                  GS.infos.underexpression="spot.list.underexpression",
                  GS.infos.samples="spot.list.samples",
                  batch.t.g.m="t.ensID.m",
                  genes.coordinates="gene.coordinates",
                  indata.mean.level="indata.gene.mean",
                  supersom.20="secLvlSOM.20.20",
                  supersom.custom="secLvlSOM.custom",
                  unique.group.colors="groupwise.group.colors")

  prefs <- list(differences.list="pairwise.comparison.list",
                dim.som1="dim.1stLvlSom",
                dim.som2="dim.2ndLvlSom",
                ensembl.dataset="database.dataset",
                ensembl.rowname.ids="database.id.type",
                feature.mean.normalization="feature.centralization",
                flip.som1="rotate.SOM.portraits",
                group.spot.core="spot.coresize.groupmap",
                group.spot.threshold="spot.threshold.groupmap",
                rotate.som1="rotate.SOM.portraits",
                sample.spot.cutoff="spot.threshold.samples",
                summary.spot.core="spot.coresize.modules",
                summary.spot.threshold="spot.threshold.modules")

  for (x in names(prefs)) {
    old.env$preferences[prefs[[x]]] <- old.env$preferences[x]
  }
  env <- opossom.new(old.env$preferences)
  rm(preferences, envir=old.env)

  for (x in ls(old.env)) {
    y <- x

    if (x %in% names(globals)) {
      y <- globals[[x]]
    }
    assign(y, get(x, envir=old.env), envir=env)
  }

  env$colramp <-
    colorRampPalette(c("darkblue", "blue", "lightblue", "green",
                       "yellow", "red", "darkred"))

  save(env, file=fout)
}

# usage: migrate("Charles Stroma 18.RData", "Stroma.RData")
#        load("Stroma.RData")
#        ls(env)
