pipeline.summarySheetsModules <- function()
{
  environment(modules.CSV.sheets) <- environment()
  environment(modules.report.sheets) <- environment()
  environment(modules.profiles) <- environment()
  environment(modules.chromosomes) <- environment()
  environment(modules.relations) <- environment()
  
  dirname <- paste(files.name, "- Results/Summary Sheets - Modules")
  dir.create(dirname, showWarnings=FALSE, recursive=TRUE)

  #### Overexpression Spots ####
  
  dirname <- file.path( paste(files.name, "- Results"),"Summary Sheets - Modules","Overexpression Spots" )
  util.info("Writing:", file.path(dirname, "*.pdf"))
  dir.create(dirname, showWarnings=FALSE, recursive=TRUE)

  modules.report.sheets(spot.list=spot.list.overexpression, main="Overexpression Spots", path=file.path(dirname,"Report.pdf") )
  modules.profiles(spot.list=spot.list.overexpression, main="Overexpression Spots", path=file.path(dirname,"Profiles.pdf") )
  modules.chromosomes(spot.list=spot.list.overexpression, main="Overexpression Spots", path=file.path(dirname,"Chromosomes.pdf") )
  modules.relations(spot.list=spot.list.overexpression, main="Overexpression Spots", path=file.path(dirname,"Relations.pdf") )
  
  
  #### Underexpression Spots ####
  
  dirname <- file.path( paste(files.name, "- Results"),"Summary Sheets - Modules","Underexpression Spots" )
  util.info("Writing:", file.path(dirname, "*.pdf"))
  dir.create(dirname, showWarnings=FALSE, recursive=TRUE)
  
  modules.report.sheets(spot.list=spot.list.underexpression, main="Underexpression Spots", path=file.path(dirname,"Report.pdf") )
  
  
  #### Correlation Cluster ####
  
  dirname <- file.path( paste(files.name, "- Results"),"Summary Sheets - Modules","Correlation Cluster" )
  util.info("Writing:", file.path(dirname, "*.pdf"))
  dir.create(dirname, showWarnings=FALSE, recursive=TRUE)
  
  modules.report.sheets(spot.list=spot.list.correlation, main="Correlation Cluster", path=file.path(dirname,"Report.pdf") )
  
  
  #### K-Means Cluster ####
  
  dirname <- file.path( paste(files.name, "- Results"),"Summary Sheets - Modules","K-Means Cluster" )
  util.info("Writing:", file.path(dirname, "*.pdf"))
  dir.create(dirname, showWarnings=FALSE, recursive=TRUE)
  
  modules.report.sheets(spot.list=spot.list.kmeans, main="K-Means Cluster", path=file.path(dirname,"Report.pdf") )
  modules.profiles(spot.list=spot.list.kmeans, main="K-Means Cluster", path=file.path(dirname,"Profiles.pdf") )
  modules.chromosomes(spot.list=spot.list.kmeans, main="K-Means Cluster", path=file.path(dirname,"Chromosomes.pdf") )
  modules.relations(spot.list=spot.list.kmeans, main="K-Means Cluster", path=file.path(dirname,"Relations.pdf") )
  
  
  #### D-Clusters ####
  
  dirname <- file.path( paste(files.name, "- Results"),"Summary Sheets - Modules","D-Cluster" )
  util.info("Writing:", file.path(dirname, "*.pdf"))
  dir.create(dirname, showWarnings=FALSE, recursive=TRUE)
  
  modules.report.sheets(spot.list=spot.list.dmap, main="D-Cluster", path=file.path(dirname,"Report.pdf") )
  modules.profiles(spot.list=spot.list.dmap, main="D-Cluster", path=file.path(dirname,"Profiles.pdf") )
  modules.chromosomes(spot.list=spot.list.dmap, main="D-Cluster", path=file.path(dirname,"Chromosomes.pdf") )
  modules.relations(spot.list=spot.list.dmap, main="D-Cluster", path=file.path(dirname,"Relations.pdf") )
  
  
  #### Group Overexpression Spots ####
  
  if (length(unique(group.labels)) > 1)
  {
    dirname <- file.path( paste(files.name, "- Results"),"Summary Sheets - Modules","Group Overexpression Spots" )
    util.info("Writing:", file.path(dirname, "*.pdf"))
    dir.create(dirname, showWarnings=FALSE, recursive=TRUE)
    
    modules.report.sheets(spot.list=spot.list.group.overexpression, main="Group Overexpression Spots", path=file.path(dirname,"Report.pdf") )
    modules.profiles(spot.list=spot.list.group.overexpression, main="Group Overexpression Spots", path=file.path(dirname,"Profiles.pdf") )
    modules.chromosomes(spot.list=spot.list.group.overexpression, main="Group Overexpression Spots", path=file.path(dirname,"Chromosomes.pdf") )
    modules.relations(spot.list=spot.list.group.overexpression, main="Group Overexpression Spots", path=file.path(dirname,"Relations.pdf") )
  }

    
  #### module gene lists CSV sheets ####
   
  dirname <- file.path(output.paths["CSV"], "Spot Lists")
  util.info("Writing:", file.path(dirname, "*.csv"))
  dir.create(dirname, showWarnings=FALSE, recursive=TRUE)
  
  modules.CSV.sheets(spot.list=spot.list.overexpression, main="Overexpression Spots", path=dirname )
  modules.CSV.sheets(spot.list=spot.list.kmeans, main="K-Means Clusters", path=dirname )
  modules.CSV.sheets(spot.list=spot.list.dmap, main="D-Clusters", path=dirname )
  
  if (length(unique(group.labels)) > 1)
  {
    modules.CSV.sheets(spot.list=spot.list.group.overexpression, main="Group Overexpression Spots", path=dirname )
  }
  
}