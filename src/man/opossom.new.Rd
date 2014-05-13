\name{opossom.new}
\alias{opossom.new}
\title{
  Initiaize the oposSOM pipeline.
}
\description{
  This function initializes the oposSOM environment and sets the preferences.
}
\usage{
  opossom.new(preferences)
}
\arguments{
  \item{preferences}{
    list with the following optional members:

    \itemize{
      \item{dataset.name}: name of the dataset
      \item{error.model}: error model to use (character \{\code{"all.samples"},\code{"replicates"},\code{"groups"}\})
      \item{dim.1stLvlSom}: dimensions of primary SOM in x and y range (numerical)
      \item{dim.2ndLvlSom}: dimensions of second level SOM in x and y range (numerical)
      \item{training.extension}: factor to extend the number of iterations in SOM trainin (numerical)
      \item{rotate.SOM.portraits}: number of roations of the primary SOM in counter-clockwise fashion (numerical \{0,1,2,3\})
      \item{flip.SOM.portraits}: mirroring the primary SOM along the bottom-left to top-right diagonal (boolean)
      \item{database.dataset}: type of ensemble dataset addressed with biomaRt interface (character)
      \item{database.id.type}: type of rowname identifier in biomaRt database (character)
      \item{geneset.analysis}: enables geneset analysis (boolean)
      \item{geneset.analysis.exact}: enables p-value and fdr calculation in geneset analysis (boolean)
      \item{max.parallel.cores}: number of CPU cores used by parallelization (numerical)
      \item{spot.threshold.samples}: spot detection threshold in single sample portraits (numerical)
      \item{summary.spot.core}: spot detection in summary maps , minimum size (numerical)
      \item{summary.spot.threshold}: spot detection in summary maps, expression threshold (numerical)
      \item{spot.coresize.groupmap}: spot detection in group-specific summary maps , minimum size (numerical)
      \item{spot.threshold.groupmap}: spot detection in group-specific summary maps, expression threshold (numerical)
      \item{feature.centralization}: enables centralization of the features (boolean)
      \item{sample.quantile.normalization}: enables quantile normalization to the samples (boolean)
      \item{pairwise.comparison.list}: group list for pairwise analyses (list, see below)
    }
  }
}
\value{
  A new oposSOM environment you can modify and pass to \code{opossom.run}.
}
\author{
Author: Henry Wirth <wirth@izbi.uni-leipzig.de>

Maintainer: Martin Kalcher <mkalcher@porkbox.net>
}
\examples{
opossom <- opossom.new(list(dataset.name="Example",
                            error.model="all.samples",
                            dim.1stLvlSom=10,
                            dim.2ndLvlSom=10,
                            pairwise.comparison.list=list(
                              list("groupA"=c("sample1", "sample2"),
                                   "groupB"=c("sample3", "sample4")))))
}
