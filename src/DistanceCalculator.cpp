#include <Rcpp.h>
#include <RcppParallel.h>

using namespace Rcpp;

struct EuclideanDistanceCalculator: RcppParallel::Worker{
  const RcppParallel::RMatrix<double> inputMatrix_;
  RcppParallel::RVector<double> resultED_;

  EuclideanDistanceCalculator(const Rcpp::NumericMatrix& inputMatrix,
                              Rcpp::NumericVector& resultED)
    : inputMatrix_(inputMatrix), resultED_(resultED){
  }

  void operator()(std::size_t begin, std::size_t end){
    for (std::size_t i = begin; i < end; i++){
      RcppParallel::RMatrix<double>::Row row = inputMatrix_.row(i);

      double partialEuclideanDistance = 0;
      for (unsigned int j = 0; j < row.length(); j++){
         partialEuclideanDistance += row[j] * row[j];
      }
      resultED_[i] = partialEuclideanDistance;
    }
  }
};

// [[Rcpp::export]]
void calculateEuclideanDistances(const Rcpp::NumericMatrix& deltaMatrix, Rcpp::NumericVector& resultEuclideanDistances2){
  int length = deltaMatrix.nrow();
  EuclideanDistanceCalculator edc(deltaMatrix, resultEuclideanDistances2);
  RcppParallel::parallelFor(0, length, edc);
}

