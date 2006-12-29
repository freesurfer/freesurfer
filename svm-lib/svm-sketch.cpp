/**
 * @file  svm-sketch.cpp
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 02:09:17 $
 *    $Revision: 1.3 $
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


#include "svm-sketch.h"

using namespace std;

bool Sketch::init(const SvmRealVector& alpha, double b, const DoubleMatrix& kernelMatrix,
                  int posCount, int negCount, double threshold) {
  if ( alpha.size() != posCount+negCount) {
    cerr << "Sketch error: the number of alphas (" << alpha.size()
    << ") does not " << "match the size of the data set ("
    << posCount+negCount << ").\n";
    setValid(false);
    return false;
  }

  // If no threshold was given, find the maximal alpha and
  // use it to set the threshold. Right now, it is set to 10^-5 of
  // the largest alpha.

  int dataCount = posCount+negCount;

  if ( threshold < 0 ) {
    threshold = 0;
    for ( int i = 0; i < dataCount; i++ )
      if ( threshold < fabs(alpha[i]) )
        threshold = fabs(alpha[i]);
    threshold *= 1e-5;
  }


  // Count non-zero alphas
  int count = 0;
  for ( int i = 0; i < dataCount; i++ )
    if ( fabs(alpha[i]) > threshold )
      count++;

  // Init the structures
  init(count);

  // Copy alphas and indices
  int label = 1;
  for ( int i = 0, index = 0; i < dataCount; i++ ) {
    if ( i == posCount)
      label = -1;
    if ( fabs(alpha[i]) > threshold ) {
      _svIndex[index] = i;
      _alpha[index] = alpha[i]*label;
      index++;
    }
  }

  _b = b;
  _kernelMatrix = kernelMatrix;

  setValid();
  return true;
}


bool Sketch::init(const SvmRealVector& alpha, const IntVector& dataIndex, double b,
                  const DoubleMatrix& kernelMatrix, int posCount, int negCount, double threshold) {
  if ( alpha.size() != dataIndex.size() ) {
    cerr << "Sketch error: the number of alphas (" << alpha.size()
    << ") and the number of indices into the original array ("
    << dataIndex.size() << ") are not equal.\n";
    setValid(false);
    return false;
  }

  // First, construct all the arrays as if the indexing were right.
  if ( !init(alpha,b,kernelMatrix,posCount,negCount,threshold) )
    return false;

  // Re-assing the indices accoriding to the dataIndex array
  for ( int i = 0; i < svCount(); i++ )
    _svIndex[i] = dataIndex[_svIndex[i]];

  setValid();
  return true;
}




double Sketch::classify (int vectorIndex) const {
  if ( ! isValid() ) {
    cerr << "Sketch error: attempting to use uninitialized classifier "
    << "to classify an example.\n";
    return 0;
  }

  double res = -_b;
  for ( int i = 0; i < _alpha.size(); i++ )
    res += _alpha[i]*_kernelMatrix[vectorIndex][_svIndex[i]];
  return res;
}




bool Sketch::read (FILE *f, bool binary) {
  setValid(false);

  if ( !_kernel.read(f) )
    return false;

  int count;
  if ( fscanf(f,"%lf\n%d\n", &_b, &count ) != 2)
    return false;

  init(count);
  if ( !_svIndex.read(f) )
    return false;

  if ( !_alpha.read(f) )
    return false;

  if ( fscanf(f,"%d\n", &count ) != 1)
    return false;
  _kernelMatrix.init(count,count);
  if ( !_kernelMatrix.read(f) )
    return false;

  setValid();
  return true;
}




bool Sketch::write (FILE *f, bool binary) const {
  if ( !_kernel.write(f) )
    return false;

  fprintf(f,"%16.12g\n%d\n", _b, svCount());
  if ( !_svIndex.write(f) )
    return false;

  if ( !_alpha.write(f) )
    return false;

  fprintf(f,"%d\n", _kernelMatrix.rows());
  if ( !_kernelMatrix.write(f) )
    return false;

  return true;
}






