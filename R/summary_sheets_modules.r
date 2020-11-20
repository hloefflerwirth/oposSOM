pipeline.summarySheetsModules <- function(env)
{
  dirname <- paste(env$files.name, "- Results/Summary Sheets - Modules")
  dir.create(dirname, showWarnings=FALSE, recursive=TRUE)

  #### Overexpression Spots ####
  
  dirname <- file.path( paste(env$files.name, "- Results"),"Summary Sheets - Modules","Overexpression Spots" )
  util.info("Writing:", file.path(dirname, "*.pdf"))
  dir.create(dirname, showWarnings=FALSE, recursive=TRUE)

  modules.report.sheets(env=env, spot.list=env$spot.list.overexpression, main="Overexpression Spots", path=file.path(dirname,"Report.pdf") )
  modules.profiles(env=env, spot.list=env$spot.list.overexpression, main="Overexpression Spots", path=file.path(dirname,"Profiles.pdf") )
  modules.chromosomes(env=env, spot.list=env$spot.list.overexpression, main="Overexpression Spots", path=file.path(dirname,"Chromosomes.pdf") )
  modules.relations(env=env, spot.list=env$spot.list.overexpression, main="Overexpression Spots", path=file.path(dirname,"Relations.pdf") )
  
  
  #### Underexpression Spots ####
  
  dirname <- file.path( paste(env$files.name, "- Results"),"Summary Sheets - Modules","Underexpression Spots" )
  util.info("Writing:", file.path(dirname, "*.pdf"))
  dir.create(dirname, showWarnings=FALSE, recursive=TRUE)
  
  modules.report.sheets(env=env, spot.list=env$spot.list.underexpression, main="Underexpression Spots", path=file.path(dirname,"Report.pdf") )
  
  
  #### Correlation Cluster ####
  
  dirname <- file.path( paste(env$files.name, "- Results"),"Summary Sheets - Modules","Correlation Cluster" )
  util.info("Writing:", file.path(dirname, "*.pdf"))
  dir.create(dirname, showWarnings=FALSE, recursive=TRUE)
  
  modules.report.sheets(env=env, spot.list=env$spot.list.correlation, main="Correlation Cluster", path=file.path(dirname,"Report.pdf") )
  
  
  #### K-Means Cluster ####
  
  dirname <- file.path( paste(env$files.name, "- Results"),"Summary Sheets - Modules","K-Means Cluster" )
  util.info("Writing:", file.path(dirname, "*.pdf"))
  dir.create(dirname, showWarnings=FALSE, recursive=TRUE)
  
  modules.report.sheets(env=env, spot.list=env$spot.list.kmeans, main="K-Means Cluster", path=file.path(dirname,"Report.pdf") )
  modules.profiles(env=env, spot.list=env$spot.list.kmeans, main="K-Means Cluster", path=file.path(dirname,"Profiles.pdf") )
  modules.chromosomes(env=env, spot.list=env$spot.list.kmeans, main="K-Means Cluster", path=file.path(dirname,"Chromosomes.pdf") )
  modules.relations(env=env, spot.list=env$spot.list.kmeans, main="K-Means Cluster", path=file.path(dirname,"Relations.pdf") )
  
  
  #### D-Clusters ####
  
  dirname <- file.path( paste(env$files.name, "- Results"),"Summary Sheets - Modules","D-Cluster" )
  util.info("Writing:", file.path(dirname, "*.pdf"))
  dir.create(dirname, showWarnings=FALSE, recursive=TRUE)
  
  modules.report.sheets(env=env, spot.list=env$spot.list.dmap, main="D-Cluster", path=file.path(dirname,"Report.pdf") )
  modules.profiles(env=env, spot.list=env$spot.list.dmap, main="D-Cluster", path=file.path(dirname,"Profiles.pdf") )
  modules.chromosomes(env=env, spot.list=env$spot.list.dmap, main="D-Cluster", path=file.path(dirname,"Chromosomes.pdf") )
  modules.relations(env=env, spot.list=env$spot.list.dmap, main="D-Cluster", path=file.path(dirname,"Relations.pdf") )
  
  
  #### Group Overexpression Spots ####
  
  if (length(unique(env$group.labels)) > 1)
  {
    dirname <- file.path( paste(env$files.name, "- Results"),"Summary Sheets - Modules","Group Overexpression Spots" )
    util.info("Writing:", file.path(dirname, "*.pdf"))
    dir.create(dirname, showWarnings=FALSE, recursive=TRUE)
    
    modules.report.sheets(env=env, spot.list=env$spot.list.group.overexpression, main="Group Overexpression Spots", path=file.path(dirname,"Report.pdf") )
    modules.profiles(env=env, spot.list=env$spot.list.group.overexpression, main="Group Overexpression Spots", path=file.path(dirname,"Profiles.pdf") )
    modules.chromosomes(env=env, spot.list=env$spot.list.group.overexpression, main="Group Overexpression Spots", path=file.path(dirname,"Chromosomes.pdf") )
    modules.relations(env=env, spot.list=env$spot.list.group.overexpression, main="Group Overexpression Spots", path=file.path(dirname,"Relations.pdf") )
  }

    
  #### module gene lists CSV sheets ####
   
  dirname <- file.path(env$output.paths["CSV"], "Spot Lists")
  util.info("Writing:", file.path(dirname, "*.csv"))
  dir.create(dirname, showWarnings=FALSE, recursive=TRUE)
  
  modules.CSV.sheets(env=env, spot.list=env$spot.list.overexpression, main="Overexpression Spots", path=dirname )
  modules.CSV.sheets(env=env, spot.list=env$spot.list.kmeans, main="K-Means Cluster", path=dirname )
  modules.CSV.sheets(env=env, spot.list=env$spot.list.dmap, main="D-Cluster", path=dirname )
  
  if (length(unique(env$group.labels)) > 1)
  {
    modules.CSV.sheets(env=env, spot.list=env$spot.list.group.overexpression, main="Group Overexpression Spots", path=dirname )
  }
  
}