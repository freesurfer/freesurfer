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
// Name: Model
//
// This file defines a model: the support vectors and their alphas.
// A model is only valid if a kernel is specified.
//
// I/O format: the first row contains the number of the support vectors,
// the dimensionality of the data and the thresold b. This is
// followed by a row of alphas and the support vectors in a binary format
// (to reduce the space).
//
// Warning: the files are typically large, as they contain the support vectors.
//
//
//  Polina Golland polina@ai.mit.edu
//
///////////////////////////////////////////////////////////////////////////


#ifndef __SVM_MODEL_H__
#define __SVM_MODEL_H__

#include "svm-sketch.h"


class Model: public Sketch {
  SvmRealMatrix _svData;


public:

  // Constructor
  Model() : Sketch() {}


  // Copy the support vectors accoridng to the svIndices.
  bool copyData (const SvmReal * const* data, int rows, int cols);
  bool copyData (const SvmRealMatrix& data) {
    return copyData(data.data(), data.rows(), data.cols() );
  }


  // Evaluate the classification value
  double classify (const SvmRealVector& x) const {
    return classify(x.data());
  };
  double classify (const SvmReal* x) const;


  // First order derivatives. The scalar version returns 0 if the
  // model has not been initialized successfully.
  // The array version returns the status flag.
  SvmReal d10(int index, const SvmReal* x) const;
  bool d10(SvmReal* res, const SvmReal* x) const;
  bool d10(SvmRealVector& res, const SvmRealVector& x) const {
    return d10(res.data(),x.data());
  }



  // Access to model parameters

  int svDim() const {
    return _svData.cols();
  }

  const SvmRealMatrix& svData() const {
    return _svData;
  }

  const SvmReal* svData(int i) const {
    return _svData[i];
  }

  SvmReal svData(int i, int j) const {
    return _svData[i][j];
  }


  // I/O
  bool read(FILE *f, bool binary = true);
  bool write(FILE *f, bool binary = true) const;
};



#endif // __SVM_MODEL_H__































