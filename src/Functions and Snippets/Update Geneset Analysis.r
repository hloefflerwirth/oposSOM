
source("Functions/f - Import.r")


dir.create( paste( files.name, "- Results" ), showWarnings=F )
dir.create( paste( files.name, "- Results/CSV Sheets" ), showWarnings=F )
dir.create( paste( files.name, "- Results/Geneset Analysis" ), showWarnings=F )
dir.create( paste( files.name, "- Results/3rd lvl Spot Analysis" ), showWarnings=F )


geneset.custom.selection = c(  "H.Tiss", "Lymphoma", "MMML CGS", "Pathw Act",  "Cancer", "Disease",  "TF",  "Glio",  "GSEA C2", "miRNA target", "miRNA target starBase", "miRNA 3UTR", "TF Tissue", "miRNA Disease", "Toxic" )

source( "Functions and Snippets/f - Load Custom Genesets.r" )
preferences$geneset.custom.list		= geneset.custom.list

preferences$geneset.analysis = T
preferences$geneset.analysis.exact = T



source("Functions/f - Prepare Annotation.r")

source("Functions/f - Geneset Statistic Samples.r")
source("Functions/f - Geneset Statistic Integral.r")
source("Functions/f - Geneset Overviews.r")
source("Functions/f - Geneset Profiles + Maps.r")

source("Functions/f - Gene Lists.r")

source("Functions/f - Summary Sheets Samples.r")
source("Functions/f - Summary Sheets Integral.r")

source("Functions/f - 3rd lvl Summary Sheets.r")	

source("Functions/f - HTML Geneset Analysis.r")

source("Functions/f - Workspace Cleanup.r")
save.image( paste( files.name, ".RData" , sep="" ) )

source("Functions/f - Group Analyses.r")
source("Functions/f - Difference Analyses.r")


