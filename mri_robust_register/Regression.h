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
 *    $Date: 2010/01/14 19:41:04 $
 *    $Revision: 1.13 $
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
#include <cassert>
#include <vnl/vnl_vector.h>
#include <vnl/vnl_matrix.h>


class Regression
{
public:

  // constructor initializing A and B
  Regression(vnl_matrix< double > & Ap,vnl_vector< double > & bp):
      A(&Ap), b(&bp),lasterror(-1),lastweight(-1),lastzero(-1),verbose(1),floatsvd(false)
  {};
  // constructor initializing B (for simple case where x is single variable and A is (...1...)^T
  Regression(vnl_vector< double > & bp):
      A(NULL), b(&bp),lasterror(-1),lastweight(-1),lastzero(-1),verbose(1),floatsvd(false)
  {};

  // Robust solver
  vnl_vector< double > getRobustEst(double sat =  SATr, double sig =  1.4826);
  // Robust solver (returning also the sqrtweights)
  vnl_vector< double > getRobustEstW(vnl_vector< double >&w, double sat =  SATr, double sig =  1.4826);

  // Least Squares
  vnl_vector < double > getLSEst ();
  vnl_vector < double > getWeightedLSEst (const vnl_vector< double > & sqrtweights);
	vnl_vector < double > getWeightedLSEstFloat(const vnl_vector< double > & sqrtweights);


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
	void setFloatSvd(bool b)
	{ 
	   floatsvd = b;
	}

  void plotPartialSat(const std::string& fname);

protected:

  vnl_vector< double > getRobustEstWAB(vnl_vector< double >&w, double sat =  SATr, double sig =  1.4826);
  double  getRobustEstWB (vnl_vector< double >&w, double sat =  SATr, double sig =  1.4826);

  double getSigmaMAD (const vnl_vector< double >& r, double d = 1.4826);
  double VectorMedian(const vnl_vector< double >& v);

  void getSqrtTukeyDiaWeights(const vnl_vector< double >& r, vnl_vector< double > &w, double sat =  SATr );
  void getTukeyBiweight(const vnl_vector< double >& r, vnl_vector< double > &w, double sat =  SATr);
  double  getTukeyPartialSat(const vnl_vector< double >& r, double sat =  SATr);

private:
  vnl_matrix< double > * A;
  vnl_vector< double > * b;
  double lasterror,lastweight,lastzero;
	int verbose;
	bool floatsvd;
};


#endif
