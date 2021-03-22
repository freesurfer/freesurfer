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



