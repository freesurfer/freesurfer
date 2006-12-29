/**
 * @file  svm-similarity.cpp
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 02:09:17 $
 *    $Revision: 1.3 $
 *
 * Copyright (C) 2002-2007,
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


#include "svm-similarity.h"

using namespace std;

double Distance (const SvmReal* v1, const SvmReal* v2, int n) {
  double res = 0;
  for ( int i = 0; i < n; i++ )
    res += (v1[i]-v2[i])*(v1[i]-v2[i]);
  return res;
}

double Distance (const SvmRealVector& v1, const SvmRealVector& v2) {
  if ( v1.size() != v2.size() ) {
    cerr << "SVM-similarity error: vector lengths are unequal: " << v1.size()
    << " and " << v2.size() << ". Returning -1.\n";
    return -1;
  }
  return Distance(v1.data(),v2.data(),v1.size());
}



void Distance (double** dist, const SvmReal* const* data, int rows, int cols) {
  for ( int i = 0; i < rows; i++ )
    for ( int j = 0; j < i; j++ )
      dist[i][j] = dist[j][i] = Distance(data[i],data[j],cols);

  for ( int i = 0; i < rows; i++ )
    dist[i][i] = 0;
}


void Distance (DoubleMatrix& dist, const SvmRealMatrix& data) {
  dist.init(data.rows(),data.rows());
  Distance(dist.data(), data.data(), data.rows(), data.cols());
}



double Product (const SvmReal* v1, const SvmReal* v2, int n) {
  double res = 0;
  for ( int i = 0; i < n; i++ )
    res += v1[i]*v2[i];
  return res;
}


double Product (const SvmRealVector& v1, const SvmRealVector& v2) {
  if ( v1.size() != v2.size() ) {
    cerr << "SVM-similarity error: vector lengths are unequal: " << v1.size()
    << " and " << v2.size() << ". Returning -1.\n";
    return -1;
  }
  return Product(v1.data(),v2.data(),v1.size());
}



void Product (double** prod, const SvmReal* const* data, int rows, int cols) {
  for ( int i = 0; i < rows; i++ )
    for ( int j = 0; j < i; j++ )
      prod[j][i] = prod[i][j] = Product(data[i],data[j],cols);

  for ( int i = 0; i < rows; i++ )
    prod[i][i] = Product(data[i],data[i],cols);
}


void Product (DoubleMatrix& prod, const SvmRealMatrix& data) {
  prod.init(data.rows(),data.rows());
  Product(prod.data(), data.data(), data.rows(), data.cols());
}






