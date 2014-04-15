
source("lib/source/import.r")


dir.create( paste( files.name, "- Results" ), showWarnings=F )
dir.create( paste( files.name, "- Results/CSV Sheets" ), showWarnings=F )
dir.create( paste( files.name, "- Results/Geneset Analysis" ), showWarnings=F )
dir.create( paste( files.name, "- Results/3rd lvl Spot Analysis" ), showWarnings=F )


geneset.custom.selection = c(  "H.Tiss", "Lymphoma", "MMML CGS", "Pathw Act",  "Cancer", "Disease",  "TF",  "Glio",  "GSEA C2", "miRNA target", "miRNA target starBase", "miRNA 3UTR", "TF Tissue", "miRNA Disease", "Toxic" )

source("lib/source/load_custom_genesets.r")
preferences$geneset.custom.list		= geneset.custom.list

preferences$geneset.analysis = T
preferences$geneset.analysis.exact = T



source("lib/source/prepare_annotation.r")

source("lib/source/geneset_statistic_samples.r")
source("lib/source/geneset_statistic_integral.r")
source("lib/source/geneset_overviews.r")
source("lib/source/geneset_profiles_and_maps.r")

source("lib/source/gene_lists.r")

source("lib/source/summary_sheets_samples.r")
source("lib/source/summary_sheets_integral.r")

source("lib/source/3rd_lvl_summary_sheets.r")

source("lib/source/html_geneset_analysis.r")

source("lib/source/workspace_cleanup.r")
save.image( paste( files.name, ".RData" , sep="" ) )

source("lib/source/group_analyses.r")
source("lib/source/difference_analyses.r")


