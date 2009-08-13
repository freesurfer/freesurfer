/**
 * @file Regression.h
 * @brief A class to solve overconstrained system A X = B
 *
 *   it uses either least squares (standard regression)
 *   or a robust estimator (Tukey's Biweight with iterative reweighted least
 *   squares)
 *
 */

/*
 * Original Author: Martin Reuter
 * CVS Revision Info:
 *    $Author: mreuter $
 *    $Date: 2009/08/13 23:15:05 $
 *    $Revision: 1.11 $
 *
 * Copyright (C) 2008-2009
 * The General Hospital Corporation (Boston, MA).
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 * Bug reports: analysis-bugs@nmr.mgh.harvard.edu
 *
 */
//
//
// written by Martin Reuter
// Nov. 4th ,2008
//

#ifndef Regression_H
#define Regression_H

#ifdef __cplusplus
extern "C"
{
#endif
#include "matrix.h"
#ifdef __cplusplus
}
#endif

#define SATr 4.685  // this is suggested for gaussian noise

#include <utility>
#include <string>
class Regression
{
public:

  // constructor initializing A and B
  Regression(MATRIX* Ap, MATRIX* Bp):
      A(Ap), B(Bp),lasterror(-1),lastweight(-1),lastzero(-1),verbose(1)
  {};
  // constructor initializing B (for simple case where x is single variable and A is (...1...)^T
  Regression(MATRIX* Bp):
      A(NULL), B(Bp),lasterror(-1),lastweight(-1),lastzero(-1),verbose(1)
  {};

  // Robust solver
  MATRIX* getRobustEst(double sat =  SATr, double sig =  1.4826);
  // Robust solver (returning also the weights)
  std::pair < MATRIX*, MATRIX* > getRobustEstW(double sat =  SATr, double sig =  1.4826);

  // Least Squares
  MATRIX* getLSEst    (MATRIX* outX = NULL);

  double getLastError()
  {
    return lasterror;
  };
  double getLastWeightPercent()
  {
    return lastweight;
  };
  double getLastZeroWeightPercent()
  {
    return lastzero;
  };
	void setVerbose(int v)
	{
	   verbose = v;
		 if (v < 0 ) verbose = 0;
		 if (v > 2 ) verbose = 2;
  }

  void plotPartialSat(const std::string& fname);

protected:
  std::pair < MATRIX*, MATRIX* > getRobustEstWAB(double sat =  SATr, double sig =  1.4826);
  std::pair < MATRIX*, MATRIX* > getRobustEstWB(double sat =  SATr, double sig =  1.4826);

  double getSigmaMAD(MATRIX *r, double d = 1.4826);
  double VectorMedian(MATRIX *v);

  MATRIX* getSqrtTukeyDiaWeights(MATRIX * r, double sat =  SATr, MATRIX * w = NULL);
  MATRIX* getTukeyBiweight(MATRIX* r, double sat =  SATr, MATRIX* w = NULL);
  double  getTukeyPartialSat(MATRIX* r, double sat =  SATr);

private:
  MATRIX* A;
  MATRIX* B;
  double lasterror,lastweight,lastzero;
	int verbose;
};


#endif
