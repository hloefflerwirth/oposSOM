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

  e <- new.env()
  e$metadata <- do.call(cbind, by(t(metadata), group.labels, colMeans)[unique(group.labels)])

  e$indata.original <- indata.original
  colnames(e$indata.original) <- group.labels[colnames(e$indata.original)]

  e$indata <- do.call(cbind, by(t(e$indata.original),
                                colnames(e$indata.original),
                                colMeans)[unique(group.labels)])

  e$indata.mean.level <- rowMeans(e$indata)

  if (preferences$feature.mean.normalization)
  {
    e$indata <- e$indata - e$indata.mean.level
  }

  e$group.colors <- group.colors[match(colnames(e$indata), group.labels)]
  e$group.labels <- group.labels[match(colnames(e$indata), group.labels)]
  names(e$group.labels) <- e$group.labels
  names(e$group.colors) <- e$group.labels

  e$output.paths <- c("LPE" = "",
                      "CSV" = paste(files.name, "- Results/Summary Sheets - Groups/CSV Sheets"),
                      "Summary Sheets Samples"= paste(files.name, "- Results/Summary Sheets - Groups/Reports"))

  e$preferences <- preferences
  e$preferences$error.model <- "replicates"

  capture.output({
    environment(pipeline.calcStatistics) <- e
    pipeline.calcStatistics()

    environment(pipeline.detectSpotsSamples) <- e
    pipeline.detectSpotsSamples()

    if (preferences$geneset.analysis)
    {
      environment(pipeline.genesetStatisticSamples) <- e
      pipeline.genesetStatisticSamples()
    }

    environment(pipeline.geneLists) <- e
    pipeline.geneLists()

    environment(pipeline.summarySheetsSamples) <- e
    pipeline.summarySheetsGroups()
  })
}
