//
// Regression is a class to solve overconstrained system A X = B
//   it uses either least squares (standard regression)
//   or a robust estimator (Tukey's Biweight with iterative reweighted least squares)
//
// written by Martin Reuter
// Nov. 4th ,2008
//

#ifndef Regression_H
#define Regression_H

#ifdef __cplusplus
extern "C" {
#endif
#include "matrix.h"
#ifdef __cplusplus
}
#endif

#define SAT 4.685  // this is suggested for gaussian noise

#include <utility>
class Regression 
{
  public:
  
     // constructor initializing A and B 
     Regression(MATRIX* Ap, MATRIX* Bp):
           A(Ap), B(Bp),lasterror(-1),lastweight(-1),lastzero(-1) {};

     // Robust solver
     MATRIX* getRobustEst(double sat =  SAT, double sig =  1.4826);
     // Robust solver (returning also the weights)
     std::pair < MATRIX*, MATRIX* > getRobustEstW(double sat =  SAT, double sig =  1.4826);
     
     // Least Squares
     MATRIX* getLSEst    (MATRIX* outX = NULL);
     
     double getLastError() {return lasterror;};
     double getLastWeightPercent() {return lastweight;};
     double getLastZeroWeightPercent() {return lastzero;};

  protected:
  
     double getSigmaMAD(MATRIX *r, double d = 1.4826);
     double VectorMedian(MATRIX *v);
     double kth_smallest(double a[], int n, int k);
     double quick_select(double arr[], int n, int k);

     MATRIX* getSqrtTukeyDiaWeights(MATRIX * r, double sat =  SAT, MATRIX * w = NULL);
     MATRIX* getTukeyBiweight(MATRIX* r, double sat =  SAT, MATRIX* w = NULL);
     
  private:
     MATRIX* A;
     MATRIX* B;
     double lasterror,lastweight,lastzero;
};


#endif
