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
// Name: vector-types
//
// This file defines different templae instanciations for vector and matrix.
//
//  Polina Golland polina@ai.mit.edu
//
///////////////////////////////////////////////////////////////////////////

#ifndef __SVM_VECTOR_TYPES_H__
#define __SVM_VECTOR_TYPES_H__

#include "svm-element-type.h"
#include "svm-vector.h"
#include "svm-matrix.h"

using IntVector = Vector<int>;
using DoubleVector = Vector<double>;

using IntMatrix = Matrix<int>;
using DoubleMatrix = Matrix<double>;

using SvmRealVector = Vector<SvmReal>;
using SvmRealMatrix = Matrix<SvmReal>;

#endif //  __SVM_VECTOR_TYPES_H__
