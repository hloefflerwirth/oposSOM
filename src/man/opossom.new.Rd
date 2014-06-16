\name{opossom.new}
\alias{opossom.new}
\title{
  Initialize the oposSOM pipeline.
}
\description{
  This function initializes the oposSOM environment and sets the preferences.
}
\usage{
  opossom.new(preferences)
}
\arguments{
  \item{preferences}{
    list with the following optional values:

    \itemize{
      \item{dataset.name}: name of the dataset; used for results folder (character)
      \item{dim.1stLvlSom}: dimension of primary SOM (integer, >5)
      \item{dim.2ndLvlSom}: dimensions of second level SOM (integer, >5)
      \item{training.extension}: factor to extend the number of iterations in SOM training (numerical, >0)
      \item{rotate.SOM.portraits}: number of roations of the primary SOM in counter-clockwise fashion (integer \{0,1,2,3\})
      \item{flip.SOM.portraits}: mirroring the primary SOM along the bottom-left to top-right diagonal (boolean)
      \item{database.dataset}: type of ensemble dataset addressed with biomaRt interface; use "auto" to detect parameter automatically (character)
      \item{database.id.type}: type of rowname identifier in biomaRt database; obsolete if \code{database.dataset="auto"} (character)
      \item{geneset.analysis}: perform geneset analysis (boolean)
      \item{geneset.analysis.exact}: enables p-value and fdr calculation in geneset analysis; obsolete if \code{geneset.analysis=F} (boolean)
      \item{max.parallel.cores}: number of CPU cores used by parallelization (numerical, >0)
      \item{spot.threshold.samples}: spot detection threshold in single sample portraits (numerical, between 0 and 1)
      \item{spot.coresize.modules}: spot detection in summary maps , minimum size (numerical, >0)
      \item{spot.threshold.modules}: spot detection in summary maps, expression threshold (numerical, between 0 and 1)
      \item{spot.coresize.groupmap}: spot detection in group-specific summary maps , minimum size (numerical, >0)
      \item{spot.threshold.groupmap}: spot detection in group-specific summary maps, expression threshold (numerical, between 0 and 1)
      \item{feature.centralization}: enables centralization of the features (boolean)
      \item{sample.quantile.normalization}: enables quantile normalization of the samples (boolean)
      \item{pairwise.comparison.list}: group list for pairwise analyses (list of group lists, see example)
    }
  }
}
\value{
  A new oposSOM environment which is passed to \code{opossom.run}.
}
\author{
Author: Henry Wirth <wirth@izbi.uni-leipzig.de>

Maintainer: Martin Kalcher <mkalcher@porkbox.net>
}
\examples{
env <- opossom.new(list(dataset.name="Example",
                        dim.1stLvlSom=10,
                        dim.2ndLvlSom=10,
                        training.extension=1,
                        rotate.SOM.portraits=0,
                        flip.SOM.portraits=FALSE,
                        database.dataset="auto",
                        geneset.analysis=TRUE,
                        geneset.analysis.exact=TRUE,
                        max.parallel.cores=2,
                        spot.threshold.samples=0.9,
                        spot.coresize.modules=4,
                        spot.threshold.modules=0.9,
                        spot.coresize.groupmap=4,
                        spot.threshold.groupmap=0.7,
                        feature.centralization=TRUE,
                        sample.quantile.normalization=TRUE,
                        pairwise.comparison.list=list(
                          list("groupA"=c("sample1", "sample2"),
                               "groupB"=c("sample3", "sample4")))))
}
