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


////SVM-LIB////////////////////////////////////////////////////////////////
//
// Name: SvmSimilary
//
// Dot-product and distance functions.
//
//  Polina Golland polina@ai.mit.edu
//
///////////////////////////////////////////////////////////////////////////

#ifndef __SVM_SIMILARITY_H__
#define __SVM_SIMILARITY_H__

#include "svm-vector-types.h"

double Distance (const SvmReal* v1, const SvmReal* v2, int n);
double Distance (const SvmRealVector& v1, const SvmRealVector& v2);

void Distance (double** dist, const SvmReal* const* data, int rows, int cols);
void Distance (DoubleMatrix& dist, const SvmRealMatrix& data);


double Product (const SvmReal* v1, const SvmReal* v2, int n);
double Product (const SvmRealVector& v1, const SvmRealVector& v2);

void Product (double** prod, const SvmReal* const* data, int rows, int cols);
void Product (DoubleMatrix& prod, const SvmRealMatrix& data);


#endif // __SVM_SIMILARITY_H__



