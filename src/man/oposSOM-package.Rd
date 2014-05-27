\name{oposSOM-package}
\alias{oposSOM-package}
\alias{oposSOM}
\docType{package}
\title{
Package for comprehensive analysis of transciptome data
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
Version: \tab 1.0.0\cr
Date: \tab 2014-27-05\cr
License: \tab GPL (>= 2)\cr
}
}
\author{
Author: Henry Wirth <wirth@izbi.uni-leipzig.de>

Maintainer: Martin Kalcher <mkalcher@porkbox.net>
}
\references{
Wirth, Loeffler, v.Bergen, Binder: Expression cartography of human tissues %
using self organizing maps. (BMC Bioinformatics 2011)

Wirth, v.Bergen, Binder: MALDI-typing of infectious algae of the genus %
Prototheca using SOM portraits. (Journal of microbiological methods 2012)

Wirth, v.Bergen, Binder: Mining SOM expression portraits: feature selection and %
integrating concepts of molecular function. (BioData Mining 2012)
}
\keyword{ package }
\examples{
library("oposSOM")
# Example with artificial data
pipeline <- opossom.new(list(dataset.name="Example",
                             error.model="all.samples",
                             dim.1stLvlSom=20,
                             dim.2ndLvlSom=20,
                             geneset.analysis=TRUE))

pipeline$indata <- matrix(runif(10000), 1000, 10)

pipeline$group.labels <- c(rep(c("ABC", "DEF", "GHI"), each=3), "JKL")

opossom.run(pipeline)

# Real Example - This will take several minutes
#pipeline <- opossom.new(list(dataset.name="Tissues",
#                             error.model="all.samples",
#                             dim.1stLvlSom=30,
#                             dim.2ndLvlSom=30,
#                             geneset.analysis=TRUE,
#                             spot.coresize.groupmap=5,
#                             spot.threshold.groupmap=0.75,
#                             pairwise.comparison.list=list(
#                               list("Homeostasis"=c(1, 2), "Imune System"=c(9, 10)),
#                               list("Homeostasis"=c(1, 2), "Muscle"=c(8))
#                             )))
#
#data(opossom.tissues)
#pipeline$indata <- opossom.tissues
#
#pipeline$group.labels <- c(rep("Homeostasis", 2),
#                           "Endocrine",
#                           "Digestion",
#                           "Exocrine",
#                           "Epithelium",
#                           "Reproduction",
#                           "Muscle",
#                           rep("Imune System", 2),
#                           rep("Nervous System", 2))
#
#opossom.run(pipeline)
}
