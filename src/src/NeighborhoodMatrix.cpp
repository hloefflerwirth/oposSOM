#include <Rcpp.h>
#include <RcppParallel.h>

using namespace Rcpp;

// #include <cmath>


inline int xyToRowNumber(int x, int y, int xdim){
  return x + y*xdim;
}
struct NeighborhoodMatrixCalculator : RcppParallel::Worker{
  const int xWinner_;
  const int yWinner_;
  const int somSize_;

  const unsigned int largerLimit_;
  const unsigned int smallerLimit_;

  const double radius_;
  RcppParallel::RVector<double> result_;

  NeighborhoodMatrixCalculator(const int& xWinner, const int& yWinner, const int& somSize,
                               const int& largerLimit, const int& smallerLimit, const double& radius, Rcpp::NumericVector& result)
      : xWinner_(xWinner), yWinner_(yWinner), somSize_(somSize), largerLimit_(largerLimit),
        smallerLimit_(smallerLimit), radius_(radius), result_(result){

  }

  void operator()(std::size_t begin, std::size_t end){
    int tempX = 0;
    int tempY = 0;
    double neighborhoodDistance = 0;

    for (std::size_t i = begin; i < end; i++){ //ich muss nicht durch somSize laufen, nur entlang des Quadranten
      for (std::size_t j = i; j <= largerLimit_; j++){ //das gr??ere Limit innen
        //if(i>4){std::cout<< i << " " << j << std::endl;}
        neighborhoodDistance = std::exp( (i*i + j*j) / (-2 * radius_*radius_));

        // Symmetrien ausnutzen (insgesamt 8...)
        tempX = xWinner_ + i;
        tempY = yWinner_ + j;

        if (tempX < somSize_ && tempY < somSize_){
          result_[xyToRowNumber(tempX, tempY, somSize_)] = neighborhoodDistance;
        }
        {
          tempX = xWinner_ + j;
          tempY = yWinner_ + i;
          if (tempX < somSize_ && tempY < somSize_){
            result_[xyToRowNumber(tempX, tempY, somSize_)] = neighborhoodDistance;
          }

          tempX = xWinner_ - i;
          tempY = yWinner_ - j;
          if (tempX >= 0 && tempY >= 0){
            result_[xyToRowNumber(tempX, tempY, somSize_)] = neighborhoodDistance;
          }

          tempX = xWinner_ - j;
          tempY = yWinner_ - i;
          if (tempX >= 0 && tempY >= 0){
            result_[xyToRowNumber(tempX, tempY, somSize_)] = neighborhoodDistance;
          }

          tempX = xWinner_ - i;
          tempY = yWinner_ + j;
          if (tempX >= 0 && tempY < somSize_){
            result_[xyToRowNumber(tempX, tempY, somSize_)] = neighborhoodDistance;
          }

          tempX = xWinner_ - j;
          tempY = yWinner_ + i;
          if (tempX >= 0 && tempY < somSize_){
            result_[xyToRowNumber(tempX, tempY, somSize_)] = neighborhoodDistance;
          }

          tempX = xWinner_ + i;
          tempY = yWinner_ - j;
          if (tempX < somSize_ && tempY >= 0){
            result_[xyToRowNumber(tempX, tempY, somSize_)] = neighborhoodDistance;
          }

          tempX = xWinner_ + j;
          tempY = yWinner_ - i;
          if (tempX < somSize_ && tempY >= 0){
            result_[xyToRowNumber(tempX, tempY, somSize_)] = neighborhoodDistance;
          }
        }
      }
    }
  }
};

// [[Rcpp::export]]
void calculateNeighborhoodMatrix(const int& winnerNeuronR, const int& somSize, const double& radius, Rcpp::NumericVector& resultVector){
  //Neighborhood-Matrix als Tabelle (gibt (somSize x somSize)-Matrix zur?ck)
  //Rcpp::NumericVector result(somSize*somSize);

  int winnerNeuron = winnerNeuronR - 1;
  int yWinner = (int) winnerNeuron / somSize;
  int xWinner = winnerNeuron % somSize;
  int xLimit = 0;
  int yLimit = 0;
  int largestQuadrant = 0; //zweistellige Zahl xy; 11: oben rechts, 01: oben links, 10: unten rechts; 00: unten links

  if (yWinner > somSize/2){
    largestQuadrant = 1;
  }
  if (xWinner < somSize/2){
    largestQuadrant += 10;
  }


  switch(largestQuadrant){ //am meisten Platz ist...
  case(11): //oben rechts
    xLimit = somSize - xWinner;
    yLimit = yWinner;
    break;
  case(10): //unten rechts
    xLimit = somSize - xWinner;
    yLimit = somSize - yWinner;
    break;
  case(01): //oben links
    xLimit = xWinner;
    yLimit = yWinner;
    break;
  case(00): //unten links
    xLimit = xWinner;
    yLimit = somSize - yWinner;
    break;
  }

  int* largerLimit = &xLimit;
  int* smallerLimit = &yLimit;
  if (xLimit < yLimit){
    largerLimit = &yLimit;
    smallerLimit = &xLimit;
  }

  NeighborhoodMatrixCalculator nmc(xWinner, yWinner, somSize, *largerLimit, *smallerLimit, radius, resultVector);
  parallelFor(0, *smallerLimit+1, nmc);
}



struct MatrixToCodebookMatrixConverter : RcppParallel::Worker{

  const RcppParallel::RVector<double> inputMatrix_;
  //const int xDim_;
  RcppParallel::RMatrix<double> result_;

  MatrixToCodebookMatrixConverter(const Rcpp::NumericVector& inputMatrix, Rcpp::NumericMatrix& result)
    : inputMatrix_(inputMatrix), result_(result){

  }

  void operator()(std::size_t begin, std::size_t end){
    for (std::size_t i = begin; i < end; i++){

      //std::transform(result_.row(i).begin(), result_.row(i).end(),
      //               result_.row(i).begin(), [](double d){ return inputMatrix_[i]; });
      std::fill(result_.row(i).begin(), result_.row(i).end(), inputMatrix_[i]);
      // for (std::size_t j = 0; j < xDim_; j++){
      //   //result_[i,j] = inputMatrix_[i];
      //   result_.row(i)[j] = inputMatrix_[i];
      // }
    }
  }
};

// [[Rcpp::export]]
void matrixToCodebookMatrix(const Rcpp::NumericVector& matrixAsVector, Rcpp::NumericMatrix& result){
  MatrixToCodebookMatrixConverter mtcmc(matrixAsVector, result);
  parallelFor(0, matrixAsVector.length(), mtcmc);
}
