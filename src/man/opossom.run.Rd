\name{opossom.run}
\alias{opossom.run}
\title{
Execute the oposSOM pipeline.
}
\description{
This function realizes the complete pipeline functionality: single gene
expression values are culstered to metagenes using a self-organizing map. Based
on these metagenes, visualizations (e.g. expression profiles), downstreaming
sample similarity analyses (e.g. hierarchical clustering, ICA) and functional
enrichment analyses are performed.
The results are given within a separate folder and can be browsed using the
summary HTML file.
}
\usage{
opossom.run(opossom)
}
\arguments{
  \item{opossom}{
    the opossom environment created with \code{opossom.new} according to the users' preferences
  }
}
\author{
Author: Henry Wirth <wirth@izbi.uni-leipzig.de>

Maintainer: Martin Kalcher <mkalcher@porkbox.net>
}
\examples{
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
