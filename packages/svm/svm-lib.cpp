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


#include "svm-lib.h"

extern "C" {
#include "svm-pr-loqo.h"
}

using namespace std;

///////////////////////////// Utility functions for svm train and cross-validate

void svmKernel(DoubleMatrix& kernelMatrix, const DoubleMatrix& distMatrix,
               const Kernel& kernel) {
  kernelMatrix.init(distMatrix.rows(), distMatrix.cols());
  for ( int i = 0; i < kernelMatrix.rows(); i++ )
    for ( int j = 0; j < kernelMatrix.cols(); j++ )
      kernelMatrix[i][j] = kernel(distMatrix[i][j]);
}


bool checkSvmParams (const DoubleMatrix& distMatrix,
                     int posCount, int negCount)
// This function performs basic sanity checks on the parameters to
// the svm functions.
{
  if ( distMatrix.rows() != distMatrix.cols() ) {
    cerr << "SVM error: the kernel table is not square ("
    << distMatrix.rows() << "x" << distMatrix.cols() << ").\n";
    return false;
  }

  if ( distMatrix.rows() != posCount+negCount) {
    cerr << "SVM error: the size of the kernel table (" << distMatrix.rows()
    << ") does not " << "match the size of the data set ("
    << posCount+negCount << ").\n";
    return false;
  }

  return true;
}



int solveSvmOptimization (SvmRealVector& svmAlpha, double& svmB,
                          const DoubleMatrix& kernelMatrix,
                          const SvmParam& svmParam,
                          int posCount, int negCount)
// This functions sets the right structures and calls the interior
// point optimization from pr_loqo. It assumes that the size of the
// structures is correct.
{
  int n = posCount+negCount; // m = 1

  // Initialize optimization vector data
  double *c = new double[n];
  double *h = new double[n*n];
  double *a = new double[n];
  double *b = new double[1];
  double *l = new double[n];
  double *u = new double[n];

  int labelI = 1;
  for ( int i = 0; i < n; i++ ) {
    // Change the label, if it's time
    if ( i == posCount )
      labelI = -1;

    // Objective function parameters
    c[i] = -1;
    int labelJ = 1;
    for ( int j = 0; j < n; j++ ) {
      if  ( j == posCount )
        labelJ = -1;
      h[i*n+j] = labelI*labelJ*kernelMatrix[i][j];
    }

    // Equality constraint
    a[i] = labelI;

    // Inequality constraints
    l[i] = 0;
    u[i] = svmParam.C;
  }
  b[0] = 0;


  // Initialize primal and dual space
  double *primal = new double[3*n];
  double *dual = new double[1+2*n];


  // Initialize optimization parameters


  /*
   * solve the quadratic programming problem
   *
   * minimize   c' * x + 1/2 x' * H * x
   * subject to A*x = b
   *            l <= x <= u
   */


  int status = pr_loqo(n,1,c,h,a,b,l,u,primal,dual,
                       svmParam.verbose-2,svmParam.sigDig,
                       svmParam.maxIterations,
                       svmParam.optEpsilon,svmParam.C/2,0);

  delete [] c;
  delete [] h;
  delete [] a;
  delete [] b;
  delete [] l;
  delete [] u;



  // Generate the result
  bool retCode = false;
  if ( status == OPTIMAL_SOLUTION ) {
    for ( int i = 0; i < n; i++ )
      svmAlpha[i] = (SvmReal)primal[i];
    svmB = dual[0];
    retCode = true;
  }
  delete [] primal;
  delete [] dual;

  return retCode;
}




/////////////////////////// Train a classifier and save it as sketch



bool svmTrain (Sketch& sketch, const DoubleMatrix& distMatrix,
               const SvmParam& svmParam, int posCount, int negCount)
// This functions sets the right structures and calls the interior
// point optimization from pr_loqo.
{
  if ( !checkSvmParams(distMatrix,posCount,negCount) )
    return false;

  SvmRealVector alpha(distMatrix.rows());
  double b = 0;

  DoubleMatrix kernelMatrix;
  svmKernel(kernelMatrix,distMatrix,svmParam.kernel);

  if ( solveSvmOptimization(alpha,b,kernelMatrix,svmParam,posCount,negCount) ) {
    sketch.setKernel(svmParam.kernel);
    return sketch.init(alpha,b,kernelMatrix,posCount,negCount,svmParam.alphaEpsilon);
  }
  return false;
}


/////////////////////////// Cross-validate


bool svmCrossValidate (Sketch* sketchArray, const DoubleMatrix& distMatrix,
                       const SvmParam& svmParam, int posCount, int negCount)
// This functions performs leave-one-out cross-validation,
// returning the set of sketches.
{
  if ( !checkSvmParams(distMatrix,posCount,negCount) )
    return false;

  DoubleMatrix kernelMatrix;
  svmKernel(kernelMatrix,distMatrix,svmParam.kernel);

  bool retCode = true;
  int dataCount = posCount+negCount;
  int newPosCount = posCount-1, newNegCount = negCount;
  double b;
  SvmRealVector alpha(dataCount-1);
  IntVector dataIndex(dataCount-1);
  DoubleMatrix newKernelMatrix(dataCount-1,dataCount-1);
  for ( int i = 0; i < dataCount; i++ ) {
    // Change the counts if reached the second class
    if ( i == posCount ) {
      newPosCount++;
      newNegCount--;
    }

    // Create new data: missing i-th row and i-th column
    for ( int j = 0; j < i; j++ ) {
      dataIndex[j] = j;
      for ( int m = 0; m < i; m++ )
        newKernelMatrix[j][m] = kernelMatrix[j][m];
      for ( int m = i+1; m < dataCount; m++ )
        newKernelMatrix[j][m-1] = kernelMatrix[j][m];
    }
    for ( int j = i+1; j < dataCount; j++ ) {
      dataIndex[j-1] = j;
      for ( int m = 0; m < i; m++ )
        newKernelMatrix[j-1][m] = kernelMatrix[j][m];
      for ( int m = i+1; m < dataCount; m++ )
        newKernelMatrix[j-1][m-1] = kernelMatrix[j][m];
    }

    if ( i == 0 ) {
      write(dataIndex,"foo.txt");
    }

    // Train
    if ( solveSvmOptimization(alpha, b, newKernelMatrix,
                              svmParam, newPosCount, newNegCount) ) {
      sketchArray[i].setKernel(svmParam.kernel);
      retCode &= sketchArray[i].init(alpha, dataIndex, b, kernelMatrix,
                                     newPosCount, newNegCount,
                                     svmParam.alphaEpsilon);
    } else
      return false;

  }
  return retCode;
}



bool svmCrossValidate (DoubleVector& label, const DoubleMatrix& distMatrix,
                       const SvmParam& svmParam, int posCount, int negCount) {
  // Train the classifiers
  int dataCount = posCount+negCount;
  Sketch* sketchArray = new Sketch[dataCount];
  if ( !svmCrossValidate(sketchArray,distMatrix,svmParam,posCount,negCount) )
    return false;

  // Create the classification answers
  label.init(dataCount);
  for ( int i = 0; i < dataCount; i++ )
    label[i] = sketchArray[i].classify(i);

  delete [] sketchArray;
  return true;
}


/////////////////////////////    Classify

double svmClassify (const Sketch& sketch, int index) {
  return sketch.classify(index);
}



double svmClassify (const Model& model, const SvmReal* vec) {
  return model.classify(vec);
}


double svmClassify (const Model& model, const SvmRealVector& vec) {
  return model.classify(vec);
}




/////////////////////////////    Compute gradient


bool svmWeights (const Model& model, SvmReal* weights, const SvmReal* vec) {
  return model.d10(weights,vec);
}


bool svmWeights (const Model& model, SvmRealVector& weights, const SvmRealVector& vec) {
  return model.d10(weights,vec);
}











