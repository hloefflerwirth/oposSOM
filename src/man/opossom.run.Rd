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
opossom.run(env)
}
\arguments{
  \item{env}{
    the opossom environment created with \code{opossom.new} according to the users' preferences
  }
}
\examples{
# Example with artificial data
env <- opossom.new(list(dataset.name="Example",
                        dim.1stLvlSom=20))

env$indata <- matrix(runif(1000), 100, 10)

opossom.run(env)

}
