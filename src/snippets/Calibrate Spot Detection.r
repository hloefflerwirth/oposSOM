
###########################################################################

####################       overexpression spots      ########################

###########################################################################



  env$preferences$spot.threshold.modules
  env$preferences$spot.coresize.modules
  env$preferences$spot.threshold.groupmap
  env$preferences$spot.coresize.groupmap

  # preferences$summary.spot.threshold = 0.95
  # preferences$summary.spot.core = 5
  # preferences$group.spot.core = 5
  # preferences$group.spot.threshold = 0.75

  env$preferences$spot.threshold.modules = 0.9
  env$preferences$spot.coresize.modules = 2
  env$preferences$spot.threshold.groupmap = 0.88
  env$preferences$spot.coresize.groupmap = 3

  util.call(pipeline.detectSpotsIntegral,env)


#  colramp = colorRampPalette(c("darkblue","blue","lightblue","green","yellow","red","darkred"))

  par(mfrow=c(2,2), mar=c(2,2,2,2))
  image(matrix(env$spot.list.overexpression$overview.map, env$preferences$dim.1stLvlSom), col=env$colramp(1000))
  image(matrix(env$spot.list.overexpression$overview.mask, env$preferences$dim.1stLvlSom), col=env$colramp(1000))
  par(new=T)
  plot(0, type="n", axes=F, xlab="", ylab="", xlim=c(0,env$preferences$dim.1stLvlSom), ylim=c(0,env$preferences$dim.1stLvlSom), xaxs="i", yaxs="i")
  points(do.call(rbind, lapply(env$spot.list.overexpression$spots, function(x) x$position)), pch=16, cex=1, col="black")
  points(do.call(rbind, lapply(env$spot.list.overexpression$spots, function(x) x$position)), pch=1, cex=1, col="white")

  length(env$spot.list.overexpression$spots)


  image(matrix(env$spot.list.group.overexpression$overview.map, env$preferences$dim.1stLvlSom), col=env$colramp(1000))
  image(matrix(env$spot.list.group.overexpression$overview.mask, env$preferences$dim.1stLvlSom), col=env$colramp(1000))
  par(new=T)
  plot(0, type="n", axes=F, xlab="", ylab="", xlim=c(0,env$preferences$dim.1stLvlSom), ylim=c(0,env$preferences$dim.1stLvlSom), xaxs="i", yaxs="i")
  points(do.call(rbind, lapply(env$spot.list.group.overexpression$spots, function(x) x$position)), pch=16, cex=1, col="black")
  points(do.call(rbind, lapply(env$spot.list.group.overexpression$spots, function(x) x$position)), pch=1, cex=1, col="white")

  length(env$spot.list.group.overexpression$spots)







# Re-run necessary parts of the pipeline

source("R/source/import.r")


dir.create(paste(files.name, "- Results"), showWarnings=F)
dir.create(paste(files.name, "- Results/CSV Sheets"), showWarnings=F)


if (preferences$geneset.analysis)
{
  dir.create(paste(files.name, "- Results/Geneset Analysis"), showWarnings=F)
  source("R/source/geneset_statistic_integral.r")
  source("R/source/geneset_overviews.r")
  source("R/source/geneset_profiles_and_maps.r")
}

source("R/source/gene_lists.r")
source("R/source/summary_sheets_integral.r")

dir.create(paste(files.name, "- Results/3rd lvl Spot Analysis"), showWarnings=F)
source("R/source/3rd_lvl_chromosomal_enrichment.r")
source("R/source/3rd_lvl_summary_sheets.r")
source("R/source/3rd_lvl_networks.r")

source("R/source/html_integral_summary.r")

source("R/source/workspace_cleanup.r")
save.image(paste(files.name, ".RData" , sep=""))


#source("R/source/3rd_lvl_overexpression_genenet.r")
source("R/source/signature_sets.r")






