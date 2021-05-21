/**
 * @brief A class to solve overconstrained system A X = b
 *
 *   it uses either least squares (standard regression)
 *   or a robust estimator (Tukey's Biweight with iterative reweighted least
 *   squares)
 *
 */

/*
 * Original Author: Martin Reuter
 *
 * Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
 *
 * Terms and conditions for use, reproduction, distribution and contribution
 * are found in the 'FreeSurfer Software License Agreement' contained
 * in the file 'LICENSE' found in the FreeSurfer distribution, and here:
 *
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
 *
 * Reporting: freesurfer@nmr.mgh.harvard.edu
 *
 */
//
//
// written by Martin Reuter
// Nov. 4th ,2008
//
#ifndef Regression_H
#define Regression_H

#include "matrix.h"

#define SATr 4.685  // this is suggested for gaussian noise
#include <utility>
#include <string>
#include <cassert>
#include <vnl/vnl_vector.h>
#include <vnl/vnl_matrix.h>

/** \class Transform3dTranslate
 * \brief Templated class for iteratively reweighted least squares
 */
template<class T>
class Regression
{
public:

  //! Constructor initializing A and b
  Regression(vnl_matrix<T> & Ap, vnl_vector<T> & bp) :
      A(&Ap), b(&bp), lasterror(-1), lastweight(-1), lastzero(-1), verbose(1), floatsvd(false)
  {}

  //! Constructor initializing b (for simple case where x is single variable and A is (...1...)^T
  Regression(vnl_vector<T> & bp) :
      A(NULL), b(&bp), lasterror(-1), lastweight(-1), lastzero(-1), verbose(1), floatsvd(false)
  {}

  //! Robust solver
  vnl_vector<T> getRobustEst(double sat = SATr, double sig=1.4826);
  //! Robust solver (returning also the sqrtweights)
  vnl_vector<T> getRobustEstW(vnl_vector<T>&w, double sat=SATr, double sig=1.4826);

  //! Least Squares
  vnl_vector<T> getLSEst();
  //! Weighted least squares
  vnl_vector<T> getWeightedLSEst(const vnl_vector<T> & sqrtweights);
  //! Weighted least squares in float (only for the T=double version)
  vnl_vector<T> getWeightedLSEstFloat(const vnl_vector<T> & sqrtweights);

  double getLastError()
  {
    return lasterror;
  }

  double getLastWeightPercent()
  {
    return lastweight;
  }

  double getLastZeroWeightPercent()
  {
    return lastzero;
  }

  //! Set verbose level
  void setVerbose(int v)
  {
    verbose = v;
    if (v < 0)
      verbose = 0;
    if (v > 2)
      verbose = 2;
  }
  
  //! Specify if SVD is float (also in double case)
  void setFloatSvd(bool b)
  {
    floatsvd = b;
  }

  void plotPartialSat(const std::string& fname);

protected:

  vnl_vector<T> getRobustEstWAB(vnl_vector<T>&w, double sat = SATr, double sig = 1.4826);
  double getRobustEstWB(vnl_vector<T>&w, double sat = SATr, double sig = 1.4826);

  T getSigmaMAD(const vnl_vector<T>& r, T d = 1.4826);
  T VectorMedian(const vnl_vector<T>& v);

  void getSqrtTukeyDiaWeights(const vnl_vector<T>& r, vnl_vector<T> &w, double sat = SATr);
  void getTukeyBiweight(const vnl_vector<T>& r, vnl_vector<T> &w, double sat = SATr);
  double getTukeyPartialSat(const vnl_vector<T>& r, double sat = SATr);

private:
  vnl_matrix<T> * A;
  vnl_vector<T> * b;
  double lasterror, lastweight, lastzero;
  int verbose;
  bool floatsvd;
};

#include "Regression.cpp"

#endif
