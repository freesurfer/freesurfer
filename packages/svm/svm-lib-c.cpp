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


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>


#include "svm-lib.h"

extern "C" {

#include "svm-lib-c.h"


#define SVM_NULL_MODE     0
#define SVM_SKETCH_MODE   1
#define SVM_MODEL_MODE    2


  static int svmMode = SVM_NULL_MODE;
  static SvmParam svmParam;

  static DoubleMatrix svmDistMatrix;
  static Model svmModel;

  static DataParam svmDataParam;


////////////////// Utility functions

  void SVMprepareDistMatrix(SVMreal** data, int dataCount, int featureCount) {
    svmDistMatrix.init(dataCount,dataCount);
    if ( svmParam.kernel.type() == RBF_KERNEL )
      Distance(svmDistMatrix.data(),data,dataCount,featureCount);
    else
      Product(svmDistMatrix.data(),data,dataCount,featureCount);
  }




////////////////// Train using data array


  int SVMtrain(SVMreal** data, int posCount, int negCount, int featureCount) {
    svmMode = SVM_NULL_MODE;

    SVMprepareDistMatrix(data,posCount+negCount,featureCount);

    if ( !svmTrain (svmModel, svmDistMatrix, svmParam, posCount, negCount) )
      return 0;

    if ( !svmModel.copyData(data,posCount+negCount,featureCount) )
      return 0;


    svmMode = SVM_MODEL_MODE;
    return 1;
  }


////////////////// Classify using a previously initialized model

  double SVMclassify(SVMreal* vec) {
    if ( svmMode != SVM_MODEL_MODE ) {
      std::cerr << "SVMclassify error: no classifier model has been created. Must train/read in "
      << "a classifier before using this function. Zero is returned.\n";
      return 0;
    }

    return svmClassify(svmModel,vec);
  }


////////////////// Compute the weights of the model


  int SVMweights(SVMreal* weights, SVMreal* vec) {
    if ( svmMode != SVM_MODEL_MODE ) {
      std::cerr << "SVMweights error: no classifier model has been created. Must train/read in "
      << "a classifier before using this function.\n";
      return 0;
    }

    return svmWeights(svmModel,weights,vec);
  }


////////////////// Cross-validate: return the labels


  int SVMcrossValidate (double* label, SVMreal** data,
                        int posCount, int negCount, int featureCount) {
    svmMode = SVM_NULL_MODE;

    SVMprepareDistMatrix(data,posCount+negCount,featureCount);
    DoubleVector labelVec;
    if ( !svmCrossValidate(labelVec, svmDistMatrix, svmParam, posCount, negCount) )
      return 0;

    for ( int i = 0; i < labelVec.size(); i++ )
      label[i] = labelVec[i];

    svmMode = SVM_SKETCH_MODE;
    return 1;
  }



/////////////////////// Accessors


  void SVMgetParam(SVMparam* svm) {
    svm->C = svmParam.C;
    sprintf(svm->kernel,"%s",svmParam.kernel.getString());

    svm->verbose = svmParam.verbose;
    svm->alphaEpsilon =  svmParam.alphaEpsilon;
    svm->classEpsilon =  svmParam.classEpsilon;

    svm->maxIterations = svmParam.maxIterations;
    svm->sigDig = svmParam.sigDig;
    svm->optEpsilon = svmParam.optEpsilon;
  }


  void SVMsetParam(SVMparam* svm) {
    svmParam.C = svm->C;
    svmParam.kernel.parse(svm->kernel);

    svmParam.verbose = svm->verbose;
    svmParam.alphaEpsilon =  svm->alphaEpsilon;
    svmParam.classEpsilon =  svm->classEpsilon;

    svmParam.maxIterations = svm->maxIterations;
    svmParam.sigDig = svm->sigDig;
    svmParam.optEpsilon = svm->optEpsilon;
  }


  int SVMgetFeatureCount(int* count) {
    if ( svmMode != SVM_MODEL_MODE ) {
      std::cerr << "SVMgetSvCount error: no classifier model has been created. Must train/read in "
      << "a classifier before using this function.\n";
      return 0;
    }

    *count = svmModel.svDim();
    return 1;
  }



  int SVMgetSvCount(int* count) {
    if ( svmMode != SVM_MODEL_MODE ) {
      std::cerr << "SVMgetSvCount error: no classifier model has been created. Must train/read in "
      << "a classifier before using this function.\n";
      return 0;
    }

    *count = svmModel.svCount();
    return 1;
  }


  int SVMgetB(double* b) {
    if ( svmMode != SVM_MODEL_MODE ) {
      std::cerr << "SVMgetB error: no classifier model has been created. Must train/read in "
      << "a classifier before using this function.\n";
      return 0;
    }

    *b = svmModel.b();
    return 1;
  }


  int SVMgetAlphas(SVMreal* alpha) {
    if ( svmMode != SVM_MODEL_MODE ) {
      std::cerr << "SVMgetAlphas error: no classifier model has been created. Must train/read in "
      << "a classifier before using this function.\n";
      return 0;
    }

    for ( int i = 0; i < svmModel.svCount(); i++ )
      alpha[i] = svmModel.alpha(i);
    return 1;
  }


  int SVMgetSvIndex(int* svIndex) {
    if ( svmMode != SVM_MODEL_MODE ) {
      std::cerr << "SVMgetSvIndex error: no classifier model has been created. Must train/read in "
      << "a classifier before using this function.\n";
      return 0;
    }
    for ( int i = 0; i < svmModel.svCount(); i++ )
      svIndex[i] = svmModel.svIndex(i);
    return 1;
  }



  int SVMgetDistMatrix(double** distMatrix) {
    if ( svmMode == SVM_NULL_MODE ) {
      std::cerr << "SVMgetDistMatrix error: no similarity table has been created. Must "
      << "train/cross-validate or read in a classifier before using this function.\n";
      return 0;
    }

    for ( int i = 0; i < svmDistMatrix.rows(); i++ )
      for ( int j = 0; j < svmDistMatrix.rows(); j++ )
        distMatrix[i][j] = svmDistMatrix[i][j];

    return 1;
  }


/////////////////////// I/O


  int SVMreadClassifierFromFile(FILE *f, int binary) {
    if ( !svmModel.read(f, binary) )
      return 0;
    svmMode = SVM_MODEL_MODE;
    return 1;
  }


  int SVMwriteClassifierToFile(FILE *f, int binary) {
    return svmModel.write(f, binary);
  }


  int SVMreadClassifier( char* fileName, int binary) {
    if ( !read(svmModel,fileName, binary) )
      return 0;

    svmMode = SVM_MODEL_MODE;
    return 1;
  }


  int SVMwriteClassifier( char* fileName, int binary) {
    return write(svmModel,fileName, binary);
  }


  void SVMprintSvmOptions() {
    svmParam.printUsage(std::cerr);
  }

  void SVMparseSvmOptions(char* argv[], int argc) {
    svmParam.parse(argv,argc);
  }



/////////////////////////////   Help functions for data I/O


  int SVMreadData(SVMreal** data, int rows, int cols) {
    SvmRealMatrix dataMatrix(rows,cols);
    if (!svmDataParam.readData(dataMatrix) )
      return 0;
    for ( int i = 0; i < dataMatrix.rows(); i++ )
      for ( int j = 0; j < dataMatrix.cols(); j++ )
        data[i][j] = (SVMreal)dataMatrix[i][j];

    return 1;
  }


  int SVMparseDataOptions(char *argv[], int argc, int *k, int posCount, int negCount) {
    svmDataParam.init(posCount,negCount);
    return svmDataParam.parse(argv,argc,*k);
  }


  void SVMprintDataOptions() {
    DataParam::printUsage(std::cerr);
  }


  void SVMprintDataOptionHelp() {
    DataParam::printUsageHelp(std::cerr);
  }


  int SVMgetBinaryFlag() {
    return svmDataParam.binary;
  }



} // extern "C"





