/*
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






