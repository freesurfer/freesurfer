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
// Name: Sketch
//
// This file defines a sketch: a classifier based on the distance
// table and support vectors. The sketch contains the indices of the
// suport vectors instead of the vectors themselves. This allows
// for a compact representation. A sketch is only valid if a kernel
// and a distance table are specified.
//
// I/O format: the number of support vectors and the threshold b in the first
// row, followed by index-alpha pairs, one pair per row. This format makes
// it easy to read the file into matlab.
//
// Note: b is the threshold, i.e., the classifier is
//
//        y(x) = sum_i alpha_i K(x,x_i) - b.
//
//
//  Polina Golland polina@ai.mit.edu
//
///////////////////////////////////////////////////////////////////////////


#ifndef __SVM_SKETCH_H__
#define __SVM_SKETCH_H__

#include <stdio.h>
#include <iostream>

#include "svm-io-format.h"
#include "svm-vector-types.h"
#include "svm-kernel.h"


class Sketch {
  bool _valid;

  double _b;
  SvmRealVector _alpha;
  IntVector _svIndex;

  Kernel _kernel;
  DoubleMatrix _kernelMatrix;


protected:

  // Initialize the sketch to a particular size. This function is made private
  // as the callers should use more substantial versions of init below.
  void init (int count) {
    _alpha.init(count);
    _svIndex.init(count);
    _kernelMatrix.init(count,count);
  }

  // This function can only be used by tis class and its derived classes.
  // The callers can only check if the sketch is valid.
  void setValid (bool valid = true) {
    _valid = valid;
  }




public:

  Sketch() : _valid(false) {}

  // Initialize the sketch from SVM's alphas. Construct
  // _svIndex to point to the support vectors in the training data set.
  // The threshold is used to assign alpha's to zero. If no threshold
  // for alpha is given, it will be set to 1e-5 of the largest alpha.

  bool init (const SvmRealVector& alpha, double b, const DoubleMatrix& kernelMatrix,
             int posCount, int negCount, double threshold = -1);

  // Initialize the sketch from SVM's alphas. This function is used if
  // the training was done using a sub-set of the original data set.
  // dataIndex contains the indices into the original array. This is
  // useful in cross-validation, chunking, etc.

  bool init (const SvmRealVector& alpha, const IntVector& dataIndex, double b,
             const DoubleMatrix& kernelMatrix, int posCount, int negCount, double threshold = -1);


  // Classify an example from the training set. This uses the table of distance
  // values constructed at the training phase to compute the value of the
  // classifier.

  double classify (int vectorIndex) const;



  // Accessor functions
  bool isValid() const {
    return _valid;
  }

  double b() const {
    return _b;
  }

  int svCount() const {
    return _alpha.size();
  }

  const Kernel& kernel() const {
    return _kernel;
  }

  SvmReal kernel(SvmReal dist) const {
    return _kernel(dist);
  }

  double kernel (const SvmReal* v1, const SvmReal* v2, int n) const {
    return _kernel(v1,v2,n);
  }

  void setKernel(const Kernel& kernel) {
    _kernel = kernel;
  }

  const SvmRealVector& alpha() const {
    return _alpha;
  }

  SvmReal alpha(int i) const {
    return _alpha[i];
  }

  const IntVector& svIndex() const {
    return _svIndex;
  }

  int svIndex(int i) const {
    return _svIndex[i];
  }

  const DoubleMatrix& kernelMatrix() const {
    return _kernelMatrix;
  }

  double kernelMatrix(int i, int j) const {
    return _kernelMatrix[i][j];
  }


  //I/O

  bool read (FILE* f, bool binary = false);
  bool write (FILE* f, bool binary = false) const;
};


#endif // __SVM_SKETCH_H__








