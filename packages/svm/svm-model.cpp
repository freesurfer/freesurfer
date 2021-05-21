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


#include "svm-model.h"

using namespace std;


bool Model::copyData (const SvmReal* const* data, int rows, int cols) {
  if ( !isValid() )
    return false;

  _svData.init(svCount(),cols);
  for ( int i = 0; i < svCount(); i++ )
    if ( svIndex(i) < rows )
      _svData.setRow(i, data[svIndex(i)]);
    else {
      cerr << "Sketch error: the support vector index (" << svIndex(i)
      << ") is outside the data array (length of "
      << rows << ").\n";
      setValid(false);
      return false;
    }

  return true;
}



///////////////  Classify


double Model::classify (const SvmReal *x) const {
  if ( ! isValid() ) {
    cerr << "Model error: attempting to use uninitialized classifier "
    << "to classify an example.\n";
    return 0;
  }

  double res = -b();
  for ( int i = 0; i < svCount(); i++ )
    res += alpha(i)*kernel(x,svData(i),svDim());
  return res;
}


///////////////  Gradient


SvmReal Model::d10(int index, const SvmReal* x) const {
  if ( ! isValid() ) {
    cerr << "Model error: attempting to use uninitialized classifier "
    << "to compute the gradient.\n";
    return 0;
  }

  double res = 0;
  for ( int i = 0; i < svCount(); i++ )
    res += alpha(i)*kernel().d10(index,x,svData(i),svDim());
  return (SvmReal)res;
}


bool Model::d10(SvmReal* res, const SvmReal* x) const {
  if ( ! isValid() ) {
    cerr << "Model error: attempting to use uninitialized classifier "
    << "to compute the gradient.\n";
    return false;
  }

  // Initialize the result
  for ( int j = 0; j < svDim(); j++ )
    res[j] = 0;


  // y is an intermidiate result, it is a gradient of the kernel
  // evaluated at the currently considered support vector.
  SvmReal* y = new SvmReal[svDim()];
  for ( int i = 0; i < svCount(); i++ ) {
    kernel().d10(y,x,svData(i),svDim());
    for ( int j = 0; j < svDim(); j++ )
      res[j] += alpha(i)*y[j];
  }

  delete [] y;
  return true;
}



///////////////  I/O

bool Model::read (FILE* f, bool binary) {
  Sketch::read(f);

  // Use fgets and not fscanf, because fscanf can skip over a character
  // if the first character in the binary data is eol, or something like that.

  int dim;
  char header[300];
  fgets(header,300,f);
  if ( sscanf(header,"%d", &dim ) != 1 )
    return false;

  _svData.init(svCount(),dim);
  return _svData.read(f,binary);
}


bool Model::write (FILE* f, bool binary) const {
  Sketch::write(f);
  fprintf(f,"%d\n", svDim());

  return _svData.write(f,binary);
}

