////SVM-LIB////////////////////////////////////////////////////////////////
//
// Name: Vector
//
// This file defines a simple vector class: a 1D array of double elements
// in contiguous memory with binary and text i/o.
//
//  Polina Golland polina@ai.mit.edu
// 
///////////////////////////////////////////////////////////////////////////



#ifndef __SVM_ELEMENT_TYPES_H__
#define __SVM_ELEMENT_TYPES_H__

#include "svm-vector.h"
#include "svm-matrix.h"

typedef Vector<int> IntVector;
typedef Vector<double> DoubleVector;

typedef Matrix<int> IntMatrix;
typedef Matrix<double> DoubleMatrix;


typedef float SvmReal;
typedef Vector<SvmReal> SvmRealVector;
typedef Matrix<SvmReal> SvmRealMatrix;


#endif //  __SVM_ELEMENT_TYPES_H__
