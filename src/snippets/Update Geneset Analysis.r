
source("R/source/import.r")


dir.create( paste( files.name, "- Results" ), showWarnings=F )
dir.create( paste( files.name, "- Results/CSV Sheets" ), showWarnings=F )
dir.create( paste( files.name, "- Results/Geneset Analysis" ), showWarnings=F )
dir.create( paste( files.name, "- Results/3rd lvl Spot Analysis" ), showWarnings=F )


geneset.custom.selection = c(  "H.Tiss", "Lymphoma", "MMML CGS", "Pathw Act",  "Cancer", "Disease",  "TF",  "Glio",  "GSEA C2", "miRNA target", "miRNA target starBase", "miRNA 3UTR", "TF Tissue", "miRNA Disease", "Toxic" )

source("R/source/load_custom_genesets.r")
preferences$geneset.custom.list		= geneset.custom.list

preferences$geneset.analysis = T
preferences$geneset.analysis.exact = T



source("R/source/prepare_annotation.r")

source("R/source/geneset_statistic_samples.r")
source("R/source/geneset_statistic_integral.r")
source("R/source/geneset_overviews.r")
source("R/source/geneset_profiles_and_maps.r")

source("R/source/gene_lists.r")

source("R/source/summary_sheets_samples.r")
source("R/source/summary_sheets_integral.r")

source("R/source/3rd_lvl_summary_sheets.r")

source("R/source/html_geneset_analysis.r")

source("R/source/workspace_cleanup.r")
save.image( paste( files.name, ".RData" , sep="" ) )

source("R/source/group_analyses.r")
source("R/source/difference_analyses.r")


