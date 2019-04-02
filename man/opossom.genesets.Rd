\name{opossom.genesets}
\alias{opossom.genesets}
\title{Additional literature genesets}
\description{Genesets collected from publications and independent analyses.}
\usage{
data(opossom.genesets)
}
\format{
The data set is stored in RData (binary) format. Each element of the list represents one distinct gene set and contains the Ensembl-IDs of the member genes.
}
\details{
The oposSOM package allows for analysing the biological background of the samples using predefined sets of genes of known biological context. A large and diverse collection of such gene sets is automatically derived from the Gene Ontology (GO) annotation database using biomaRt interface. \code{opossom.genesets} contains more than 4,500 additional gene sets collected from Biocarta, KEGG and Reactome databases, from literature on chemical and genetic perturbations, from literature on cancer types and subtypes, and from previous analyses using the oposSOM pipeline.
}