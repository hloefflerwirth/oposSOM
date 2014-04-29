get.SD.estimate = function(data, samples, lambda=0.5)
{
  if (length(samples) > 1)
  {
    sd.g.m = apply(data[,samples], 1, sd)
    LPE.g.m = rep(NA, nrow(data))
    names(LPE.g.m) = rownames(data)

    o = order(rowMeans(data[,samples]))
    LPE.g.m[o] = Get.Running.Average(sd.g.m[o], min(200, round(nrow(data)*0.02)))

    SD2 = LPE.g.m[o]
    for (j in (length(SD2)-1):1)
    {
      if (SD2[j] < SD2[j+1]) SD2[j] = SD2[j+1]
    }
    LPE.g.m[o] = SD2

    sd.estimate = sqrt(lambda * sd.g.m^2 + (1 - lambda) * LPE.g.m^2)
  } else
  {
    sd = apply(data, 1, sd)
    LPE.g.m = rep(NA, nrow(data))
    names(LPE.g.m) = rownames(data)

    o = order(rowMeans(data))
    LPE.g.m[o] = Get.Running.Average(sd[o] , min(200, round(nrow(data)*0.02)))
    LPE.g.m[which(is.nan(LPE.g.m))] = 0.0000000001
    LPE.g.m[which(LPE.g.m == 0)] = 0.0000000001

    SD2 = LPE.g.m[o]
    for (j in (length(SD2)-1):1)
    {
      if (SD2[j] < SD2[j+1]) SD2[j] = SD2[j+1]
    }
    LPE.g.m[o] = SD2

    sd.estimate = LPE.g.m
  }

  return(sd.estimate)
}

pipeline.differenceAnalyses = function()
{
  differences.list <- list()

  if (length(unique(group.labels)) < 8)
  {
    differences.list <- apply(combn(unique(group.labels), 2), 2, function(x)
    {
      list(which(group.labels==x[1]), which(group.labels==x[2]))
    })

    names(differences.list) <-
      apply(combn(unique(group.labels), 2), 2, paste, collapse=" vs ")

    differences.list <- c(preferences$differences.list, differences.list)
  }

  if (length(differences.list) <= 0)
  {
    return()
  }

  util.info("Processing Differences Analyses")
  dir.create(paste(files.name, "- Results/Summary Sheets - Differences"), showWarnings=F)
  dir.create(paste(files.name, "- Results/Summary Sheets - Differences/CSV Sheets"), showWarnings=F)

  WAD.g.m <<- matrix(NA, nrow(indata), length(differences.list), dimnames=list(rownames(indata), names(differences.list)))
  t.g.m <<- matrix(NA, nrow(indata), length(differences.list), dimnames=list(rownames(indata), names(differences.list)))
  p.g.m <<- matrix(NA, nrow(indata), length(differences.list), dimnames=list(rownames(indata), names(differences.list)))
  fdr.g.m <<- matrix(NA, nrow(indata), length(differences.list), dimnames=list(rownames(indata), names(differences.list)))
  Fdr.g.m <<- matrix(NA, nrow(indata), length(differences.list), dimnames=list(rownames(indata), names(differences.list)))
  n.0.m <<- rep(NA, length(differences.list))
  names(n.0.m) <<- names(differences.list)
  perc.DE.m <<- rep(NA, length(differences.list))
  names(perc.DE.m) <<- names(differences.list)

  indata <<- indata
  metadata <<- metadata

  for (d in 1:length(ifferences.list))
  {
    samples.indata <-
      list(differences.list[[d]][[1]], differences.list[[d]][[2]])

    samples.indata.original <-
      list(which(colnames(indata.original) %in% colnames(indata)[samples.indata[[1]]]),
           which(colnames(indata.original) %in% colnames(indata)[samples.indata[[2]]]))

    n <- lapply(samples.indata.original, length)

    indata <<- cbind(indata, rowMeans(indata.original[,samples.indata.original[[1]],drop=F]) -
                                rowMeans(indata.original[,samples.indata.original[[2]],drop=F]))

    metadata <<- cbind(metadata, rowMeans(metadata[,samples.indata[[1]],drop=F]) -
                                    rowMeans(metadata[,samples.indata[[2]],drop=F]))

    sd.shrink.1 <- get.SD.estimate(data=indata.original, samples=samples.indata.original[[1]])
    sd.shrink.2 <- get.SD.estimate(data=indata.original, samples=samples.indata.original[[2]])


    t.g.m[,d] <<- rowMeans(indata.original[,samples.indata.original[[1]],drop=F]) -
                   rowMeans(indata.original[,samples.indata.original[[2]],drop=F])

    t.g.m[,d] <<- t.g.m[,d] / sqrt(sd.shrink.1 / n[[1]] + sd.shrink.2 / n[[2]])


    suppressWarnings({
      try.res <- try({
        fdrtool.result <- fdrtool(t.g.m[,d], verbose=F, plot=F)
      }, silent=T)
    })

    if (class(try.res) != "try-error")
    {
      p.g.m[,d] <<- fdrtool.result$pval
      fdr.g.m[,d] <<- fdrtool.result$lfdr
      Fdr.g.m[,d] <<- fdrtool.result$qval

      n.0.m[d] <<- fdrtool.result$param[1,"eta0"]
      perc.DE.m[d] <<- 1 - n.0.m[d]
    } else
    {
      p.g.m[,d] <<- order(indata[,d]) / nrow(indata)
      fdr.g.m[,d] <<- p.g.m[,d]
      Fdr.g.m[,d] <<- p.g.m[,d]

      n.0.m[d] <<- 0.5
      perc.DE.m[d] <<- 1 - n.0.m[d]
    }

    delta.e.g.m <- rowMeans(indata.original[,samples.indata.original[[1]],drop=F]) -
                     rowMeans(indata.original[,samples.indata.original[[2]],drop=F])

    w.g.m <<- (delta.e.g.m - min(delta.e.g.m)) / (max(delta.e.g.m) - min(delta.e.g.m))
    WAD.g.m[,d] <<- w.g.m * delta.e.g.m
  }

  indata <<- indata[, -c(1:length(group.labels)), drop=F]
  colnames(indata) <<- names(differences.list)

  metadata <<- metadata[, -c(1:length(group.labels)), drop=F]
  colnames(metadata) <<- names(differences.list)

  group.labels <<- names(differences.list)
  names(group.labels) <<- names(differences.list)
  group.colors <<- rep("gray20",length(differences.list))
  names(group.colors) <<- names(differences.list)

  output.paths <<- c("CSV" = paste(files.name, "- Results/Summary Sheets - Differences/CSV Sheets"),
                     "Summary Sheets Samples"= paste(files.name, "- Results/Summary Sheets - Differences/Reports"))

  environment(pipeline.detectSpotsSamples) <- environment()
  pipeline.detectSpotsSamples()

  if (preferences$geneset.analysis)
  {
    if (ncol(t.g.m) == 1)
    {
      # crack for by command, which requires >=2 columns
      t.g.m <<- cbind(t.g.m, t.g.m)
    }

    environment(pipeline.genesetStatisticSamples) <- environment()
    pipeline.genesetStatisticSamples()
  }

  environment(pipeline.geneLists) <- environment()
  pipeline.geneLists()

  environment(pipeline.summarySheetsSamples) <- environment()
  pipeline.summarySheetsSamples()
}
