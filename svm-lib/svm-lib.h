/**
 * @file  svm-lib.h
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


#ifndef __SVM_LIB_H__
#define __SVM_LIB_H__

#include "svm-io-format.h"
#include "svm-vector.h"
#include "svm-matrix.h"
#include "svm-element-type.h"
#include "svm-vector-types.h"
#include "svm-similarity.h"

#include "svm-kernel.h"
#include "svm-sketch.h"
#include "svm-model.h"
#include "svm-data-param.h"
#include "svm-param.h"




bool svmTrain (Sketch& sketch, const DoubleMatrix& distTable,
               const SvmParam& svmParam, int posCount, int negCount);


bool svmCrossValidate (DoubleVector& label, const DoubleMatrix& distTable,
                       const SvmParam& svmParam, int posCount, int negCount);

bool svmCrossValidate (Sketch* sketch, const DoubleMatrix& distTable,
                       const SvmParam& svmParam, int posCount, int negCount);



double svmClassify (const Sketch& sketch, int index);
double svmClassify (const Model& model, const SvmReal* vec);
double svmClassify (const Model& model, const SvmRealVector& vec);

bool svmWeights (const Model& model, SvmReal* weights, const SvmReal* vec);
bool svmWeights (const Model& model, SvmRealVector& weights, const SvmRealVector& vec);




#endif // __SVM_LIB_H__



