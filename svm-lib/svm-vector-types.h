/**
 * @file  svm-vector-types.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 02:09:17 $
 *    $Revision: 1.2 $
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

typedef Vector<int> IntVector;
typedef Vector<double> DoubleVector;

typedef Matrix<int> IntMatrix;
typedef Matrix<double> DoubleMatrix;


typedef Vector<SvmReal> SvmRealVector;
typedef Matrix<SvmReal> SvmRealMatrix;


#endif //  __SVM_VECTOR_TYPES_H__
