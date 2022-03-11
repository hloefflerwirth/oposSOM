\name{oposSOM-package}
\alias{oposSOM-package}
\alias{oposSOM}
\docType{package}
\title{
Comprehensive analysis of transciptome data
}
\description{
This package translates microarray expression data into metadata of reduced
dimension. It provides various sample-centered and group-centered visualizations,
sample similarity analyses and functional enrichment analyses. The underlying
SOM algorithm combines feature clustering, multidimensional scaling and dimension
reduction, along with strong visualization capabilities. It enables extraction
and description of functional expression modules inherent in the data.
The results are given within a separate folder and can be browsed using the
summary HTML file.
}
\details{
\tabular{ll}{
Package: \tab oposSOM\cr
Type: \tab Package\cr
Version: \tab 2.1.8\cr\cr
Date: \tab 2022-03-11\cr
License: \tab GPL (>= 2)\cr
}
}
\author{
Author: Henry Loeffler-Wirth <wirth@izbi.uni-leipzig.de> and Martin Kalcher <mkalcher@porkbox.net>

Maintainer: Henry Loeffler-Wirth <wirth@izbi.uni-leipzig.de>
}
\references{
Wirth, Loeffler, v.Bergen, Binder: Expression cartography of human tissues %
using self organizing maps. (BMC Bioinformatics 2011)

Wirth, v.Bergen, Binder: Mining SOM expression portraits: feature selection and %
integrating concepts of molecular function. (BioData Mining 2012)

Loeffler-Wirth, Kalcher, Binder: oposSOM: R-package for high-dimensional portraying %
of genome-wide expression landscapes on Bioconductor. (Bioinformatics 2015)
}
\keyword{ package }
\examples{
# Example with artificial data
env <- opossom.new(list(dataset.name="Example",
                        dim.1stLvlSom=20))

env$indata <- matrix(rnorm(10000), 1000, 10)

env$group.labels <- "auto"

opossom.run(env)

# Real Example - This will take several minutes
#env <- opossom.new(list(dataset.name="Tissues",
#                        dim.1stLvlSom=30,
#                        geneset.analysis=TRUE,
#                        pairwise.comparison.list=list(
#                          list("Homeostasis"=c(1, 2), "Imune System"=c(9, 10)),
#                          list("Homeostasis"=c(1, 2), "Muscle"=c(8))
#                        )))
#
#data(opossom.tissues)
#env$indata <- opossom.tissues
#
#env$group.labels <- c(rep("Homeostasis", 2),
#                      "Endocrine",
#                      "Digestion",
#                      "Exocrine",
#                      "Epithelium",
#                      "Reproduction",
#                      "Muscle",
#                      rep("Imune System", 2),
#                      rep("Nervous System", 2))
#
#opossom.run(env)
}
