/**
 * @file  svm-similarity.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:40 $
 *    $Revision: 1.3 $
 *
 * Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
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



