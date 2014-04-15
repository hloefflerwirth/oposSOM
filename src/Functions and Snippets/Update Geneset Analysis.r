
source("lib/f - Import.r")


dir.create( paste( files.name, "- Results" ), showWarnings=F )
dir.create( paste( files.name, "- Results/CSV Sheets" ), showWarnings=F )
dir.create( paste( files.name, "- Results/Geneset Analysis" ), showWarnings=F )
dir.create( paste( files.name, "- Results/3rd lvl Spot Analysis" ), showWarnings=F )


geneset.custom.selection = c(  "H.Tiss", "Lymphoma", "MMML CGS", "Pathw Act",  "Cancer", "Disease",  "TF",  "Glio",  "GSEA C2", "miRNA target", "miRNA target starBase", "miRNA 3UTR", "TF Tissue", "miRNA Disease", "Toxic" )

source( "Functions and Snippets/f - Load Custom Genesets.r" )
preferences$geneset.custom.list		= geneset.custom.list

preferences$geneset.analysis = T
preferences$geneset.analysis.exact = T



source("lib/f - Prepare Annotation.r")

source("lib/f - Geneset Statistic Samples.r")
source("lib/f - Geneset Statistic Integral.r")
source("lib/f - Geneset Overviews.r")
source("lib/f - Geneset Profiles + Maps.r")

source("lib/f - Gene Lists.r")

source("lib/f - Summary Sheets Samples.r")
source("lib/f - Summary Sheets Integral.r")

source("lib/f - 3rd lvl Summary Sheets.r")	

source("lib/f - HTML Geneset Analysis.r")

source("lib/f - Workspace Cleanup.r")
save.image( paste( files.name, ".RData" , sep="" ) )

source("lib/f - Group Analyses.r")
source("lib/f - Difference Analyses.r")


