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






