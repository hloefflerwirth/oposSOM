#include <Rcpp.h>
#include <RcppParallel.h>

using namespace Rcpp;



double softMinus(const double& x, const double& y){
  double result = 0;
  if (!Rcpp::traits::is_na<REALSXP>(x)){
    result = x - y;
  }
  return result;
}

	
struct DeltaMatrixCalculator: RcppParallel::Worker{
  const RcppParallel::RMatrix<double> inputMatrix_;
  const RcppParallel::RVector<double> inputVector_;
  RcppParallel::RMatrix<double> resultDelta_;

  DeltaMatrixCalculator(const Rcpp::NumericMatrix& inputMatrix,
                        const Rcpp::NumericVector& inputVector,
                        Rcpp::NumericMatrix& resultDelta)
    : inputMatrix_(inputMatrix), inputVector_(inputVector), resultDelta_(resultDelta){

  }

  void operator()(std::size_t begin, std::size_t end){
    for (std::size_t i = begin; i < end; i++){
      RcppParallel::RMatrix<double>::Row row = inputMatrix_.row(i);
      std::transform(
        inputVector_.begin(), inputVector_.end(),
        row.begin(),
        resultDelta_.row(i).begin(),
        softMinus
      );

    }
  }
};

struct DeltaMatrixCalculatorNoNA: RcppParallel::Worker{
  const RcppParallel::RMatrix<double> inputMatrix_;
  const RcppParallel::RVector<double> inputVector_;
  RcppParallel::RMatrix<double> resultDelta_;

  DeltaMatrixCalculatorNoNA(const Rcpp::NumericMatrix& inputMatrix,
                        const Rcpp::NumericVector& inputVector,
                        Rcpp::NumericMatrix& resultDelta)
    : inputMatrix_(inputMatrix), inputVector_(inputVector), resultDelta_(resultDelta){

  }

  void operator()(std::size_t begin, std::size_t end){
    for (std::size_t i = begin; i < end; i++){
      RcppParallel::RMatrix<double>::Row row = inputMatrix_.row(i);
      std::transform(
        inputVector_.begin(), inputVector_.end(),
        row.begin(),
        resultDelta_.row(i).begin(),
        std::minus<double>()
      );
    }
  }
};



// [[Rcpp::export]]
void calculateDelta(const Rcpp::NumericMatrix& inputMatrix, const Rcpp::NumericVector& inputVector, const bool naExist, Rcpp::NumericMatrix& resultDelta){
  if (naExist){
    DeltaMatrixCalculator dmc(inputMatrix, inputVector, resultDelta);
    RcppParallel::parallelFor(0, inputMatrix.nrow(), dmc);
  }else{
    DeltaMatrixCalculatorNoNA dmc(inputMatrix, inputVector, resultDelta);
    RcppParallel::parallelFor(0, inputMatrix.nrow(), dmc);
	}
}

