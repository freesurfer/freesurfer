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
