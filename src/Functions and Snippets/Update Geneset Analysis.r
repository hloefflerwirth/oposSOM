
source("lib/source/Import.r")


dir.create( paste( files.name, "- Results" ), showWarnings=F )
dir.create( paste( files.name, "- Results/CSV Sheets" ), showWarnings=F )
dir.create( paste( files.name, "- Results/Geneset Analysis" ), showWarnings=F )
dir.create( paste( files.name, "- Results/3rd lvl Spot Analysis" ), showWarnings=F )


geneset.custom.selection = c(  "H.Tiss", "Lymphoma", "MMML CGS", "Pathw Act",  "Cancer", "Disease",  "TF",  "Glio",  "GSEA C2", "miRNA target", "miRNA target starBase", "miRNA 3UTR", "TF Tissue", "miRNA Disease", "Toxic" )

source("Functions and Snippets/f - Load Custom Genesets.r")
preferences$geneset.custom.list		= geneset.custom.list

preferences$geneset.analysis = T
preferences$geneset.analysis.exact = T



source("lib/source/Prepare Annotation.r")

source("lib/source/Geneset Statistic Samples.r")
source("lib/source/Geneset Statistic Integral.r")
source("lib/source/Geneset Overviews.r")
source("lib/source/Geneset Profiles + Maps.r")

source("lib/source/Gene Lists.r")

source("lib/source/Summary Sheets Samples.r")
source("lib/source/Summary Sheets Integral.r")

source("lib/source/3rd lvl Summary Sheets.r")

source("lib/source/HTML Geneset Analysis.r")

source("lib/source/Workspace Cleanup.r")
save.image( paste( files.name, ".RData" , sep="" ) )

source("lib/source/Group Analyses.r")
source("lib/source/Difference Analyses.r")


