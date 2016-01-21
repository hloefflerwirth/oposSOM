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
      \item{indata}: input data matrix containing the expression values or an \code{Biobase::ExpressionSet} object (see 'Details' and 'Examples')
      \item{group.labels}: sample assignment to a distinct group, subtype or class (character; "auto" or one label for each sample; may be given with \code{indata ExpressionSet})
      \item{group.colors}: colors of the samples for diverse visualizations (character; one color for each sample; may be given with \code{indata ExpressionSet})
      \item{dim.1stLvlSom}: dimension of primary SOM; use "auto" to apply automatic size estimation (integer, >5)
      \item{dim.2ndLvlSom}: dimensions of second level SOM (integer, >5)
      \item{training.extension}: factor to extend the number of iterations in SOM training (numerical, >0)
      \item{rotate.SOM.portraits}: number of roations of the primary SOM in counter-clockwise fashion (integer \{0,1,2,3\})
      \item{flip.SOM.portraits}: mirroring the primary SOM along the bottom-left to top-right diagonal (boolean)
      \item{database.dataset}: type of ensemble dataset addressed with biomaRt interface; use "auto" to detect parameter automatically (character)
      \item{database.id.type}: type of rowname identifier in biomaRt database; obsolete if \code{database.dataset="auto"} (character)
      \item{geneset.analysis}: perform geneset analysis (boolean)
      \item{geneset.analysis.exact}: enables p-value and fdr calculation in geneset analysis; obsolete if \code{geneset.analysis=F} (boolean)
      \item{max.parallel.cores}: number of CPU cores used by parallelization (numerical, >0)
      \item{spot.coresize.modules}: spot detection in summary maps , minimum size (numerical, >0)
      \item{spot.threshold.modules}: spot detection in summary maps, expression threshold (numerical, between 0 and 1)
      \item{spot.coresize.groupmap}: spot detection in group-specific summary maps , minimum size (numerical, >0)
      \item{spot.threshold.groupmap}: spot detection in group-specific summary maps, expression threshold (numerical, between 0 and 1)
      \item{feature.centralization}: enables centralization of the features (boolean)
      \item{sample.quantile.normalization}: enables quantile normalization of the samples (boolean)
      \item{pairwise.comparison.list}: group list for pairwise analyses (list of group lists, see 'Examples') or NULL otherwise
    }
  }
}
\details{
The package then accepts the \code{indata} parameter in two formats:<br>
Firstly a simple two-dimensional numerical matrix, where the columns and rows represent the samples and genes, respectively. The expression values are usually obtained by calibration and summarization algorithms (e.g. MAS5, VSN or RMA), and transformed into logarithmic scale prior to utilizing them in the pipeline. Secondly the input data can also be given as \code{Biobase::ExpressionSet} object.
Please check the vignette for more details on the parameters.
}
\value{
  A new oposSOM environment which is passed to \code{opossom.run}.
}
\author{
Author: Henry Loeffler-Wirth <wirth@izbi.uni-leipzig.de>

Maintainer: Martin Kalcher <mkalcher@porkbox.net>
}
\examples{
env <- opossom.new(list(dataset.name="Example",
                        dim.1stLvlSom="auto",
                        dim.2ndLvlSom=10,
                        training.extension=1,
                        rotate.SOM.portraits=0,
                        flip.SOM.portraits=FALSE,
                        database.dataset="auto",
                        geneset.analysis=TRUE,
                        geneset.analysis.exact=TRUE,
                        max.parallel.cores=2,
                        spot.coresize.modules=4,
                        spot.threshold.modules=0.9,
                        spot.coresize.groupmap=4,
                        spot.threshold.groupmap=0.7,
                        feature.centralization=TRUE,
                        sample.quantile.normalization=TRUE,
                        pairwise.comparison.list=list(
                          list("groupA"=c("sample1", "sample2"),
                               "groupB"=c("sample3", "sample4")))))


# definition of indata, group.labels and group.colors
env$indata = matrix( runif(1000), 100, 10 )
env$group.labels = c( rep("class 1", 5), rep("class 2", 4), "class 3" )
env$group.colors = c( rep("red", 5), rep("blue", 4), "green" )

# alternative definition of indata, group.labels and group.colors using Biobase::ExpressionSet
library(Biobase)

env$indata = ExpressionSet( assayData=matrix(runif(1000), 100, 10),
                            phenoData=AnnotatedDataFrame(data.frame( 
                                group.labels = c( rep("class 1", 5), rep("class 2", 4), "class 3" ),
                                group.colors = c( rep("red", 5), rep("blue", 4), "green" ) ))
                          )

}
