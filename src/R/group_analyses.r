pipeline.groupAnalysis <- function()
{
  if (length(unique(group.labels)) <= 1)
  {
    return()
  }

  util.info("Processing Group-centered Analyses")
  dir.create(paste(files.name, "- Results/Summary Sheets - Groups"), showWarnings=F)
  dir.create(paste(files.name, "- Results/Summary Sheets - Groups/CSV Sheets"), showWarnings=F)

  if (preferences$geneset.analysis)
  {
    dir.create(paste(files.name, "- Results/Summary Sheets - Groups/Geneset Analysis"), showWarnings=F)

    environment(pipeline.groupSpecificGenesets) <- environment()
    pipeline.groupSpecificGenesets()
  }

  environment(pipeline.summarySheetsGroups) <- environment()
  pipeline.summarySheetsGroups()

  metadata <<- do.call(cbind, by(t(metadata), group.labels, colMeans)[unique(group.labels)])

  colnames(indata.original) <<- group.labels[colnames(indata.original)]

  indata <<- do.call(cbind, by(t(indata.original),
                               colnames(indata.original),
                               colMeans)[unique(group.labels)])

  indata.gene.mean <<- rowMeans(indata)

  if (preferences$feature.centralization)
  {
    indata <<- indata - indata.gene.mean
  }

  group.colors <<- group.colors[match(colnames(indata), group.labels)]
  group.labels <<- group.labels[match(colnames(indata), group.labels)]
  names(group.labels) <<- group.labels
  names(group.colors) <<- group.labels

  output.paths <<- c("LPE" = "",
                     "CSV" = paste(files.name, "- Results/Summary Sheets - Groups/CSV Sheets"),
                     "Summary Sheets Samples"= paste(files.name, "- Results/Summary Sheets - Groups/Reports"))

  preferences$error.model <<- "replicates"

  environment(pipeline.calcStatistics) <- environment()
  pipeline.calcStatistics()

  environment(pipeline.detectSpotsSamples) <- environment()
  pipeline.detectSpotsSamples()

  if (preferences$geneset.analysis)
  {
    environment(pipeline.genesetStatisticSamples) <- environment()
    pipeline.genesetStatisticSamples()
  }

  environment(pipeline.geneLists) <- environment()
  pipeline.geneLists()

  environment(pipeline.summarySheetsSamples) <- environment()
  pipeline.summarySheetsSamples()
}
