\name{opossom.run}
\alias{opossom.run}
\title{
Execute the oposSOM pipeline.
}
\description{
This function realizes the complete pipeline functionality: single gene
expression values are culstered to metagenes using a self-organizing map.
Based on these metagenes, visualizations (e.g. expression profiles) and
downstreaming analysis approaches (e.g. hierarchical clustering, ICA or
geneset overrepresentation analysis) are performed.
}
\usage{
opossom.run(opossom)
}
\arguments{
  \item{opossom}{
    the opossom environment created with \code{opossom.new} and adapted to
    your needs
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
}
