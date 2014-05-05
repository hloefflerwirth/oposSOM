\name{oposSOM-package}
\alias{oposSOM-package}
\alias{oposSOM}
\docType{package}
\title{
Package for analysis and visualization of metagenes
}
\description{
This package translates microarray expression data into metagenes and provides
various visualizations and second/third level analyses.
}
\details{
\tabular{ll}{
Package: \tab oposSOM\cr
Type: \tab Package\cr
Version: \tab 0.3.0\cr
Date: \tab 2014-05-05\cr
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
}
\keyword{ package }
\examples{
library("oposSOM")
# Fun Example
pipeline <- opossom.new(list(dataset.name="Example",
                             error.model="all.samples",
                             dim.som1=20,
                             dim.som2=20,
                             geneset.analysis=TRUE))

pipeline$indata <- matrix(runif(10000), 1000, 10)

pipeline$group.labels <- c(rep(c("ABC", "DEF", "GHI"), each=3), "JKL")

opossom.run(pipeline)

# Real Example - This will take some time
#pipeline <- opossom.new(list(dataset.name="Tissues",
#                             error.model="all.samples",
#                             dim.som1=30,
#                             dim.som2=30,
#                             geneset.analysis=TRUE,
#                             group.spot.core=5,
#                             group.spot.threshold=0.75,
#                             differences.list=list(
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
