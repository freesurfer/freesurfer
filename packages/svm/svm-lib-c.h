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


/*------------------------------------------------------------------------
 *
 * NAME: SVM
 *
 *
 * In all functions, the following inputs have fixed meaning:
 *
 * data - the original input data, one row per training example.
 *        The examples are arranged so that all positive examples preceed
 *        all negative examples.
 *
 * posCount - number of positive examples (rows) in the training data matrix.
 * negCount - number of negative examples (rows) in the training data matrix.
 * featureCount - dimensionality of the data (number of columns in the training
 *                data matrix).
 *
 * vec - test example, must be of the same dimensionality as data.
 *
 *
 *
 *  Polina Golland    polina@ai.mit.edu       02/10/2002
 *
 *
 *-------------------------------------------------------------------------*/

#ifndef __SVM_LIB_C_H__
#define __SVM_LIB_C_H__


/* Type of the feature vector elements. Should work with float and double. */
#include "svm-element-type.h"
typedef SvmReal SVMreal;


/* A set of parameters used by the svm optimization */
typedef struct {
  double C;                              /* soft margin constant */
  char kernel[300];                      /* kernel string, see kernel-param.h for options */

  int verbose;                           /* verbose mode */
  double alphaEpsilon;                   /* threshold under which alpha is considered zero */
  double classEpsilon ;                  /* epsilon around the decision boundary */

  int maxIterations;                     /* max. number of iterations in the optimization */
  double sigDig;                         /* primal and dual must agree to that many digits */
  double optEpsilon;                     /* epsilon around the constrains */

}
SVMparam;




/*---------------------  Train, classify and get the weights -------------------------------*/

int SVMtrain(SVMreal** data, int posCount, int negCount, int featureCount);

double SVMclassify(SVMreal* vec);

/* weights is the weights vector at point vec. If using linear classifier, vec could be NULL */
int SVMweights(SVMreal* weights, SVMreal* vec);


/*-------------------------------  Cross-validate --------------------------------------------*/

/* label contains the results of applying the classifiers to the corresponding hold-out point */
int SVMcrossValidate (double *label, SVMreal** data, int posCount, int negCount, int featureCount);


/*-----------------------------  Get/set optimization parameters  ----------------------------*/

void SVMgetParam(SVMparam* svm);
void SVMsetParam(SVMparam* svm);



/*-----------------------------  Get the classifier parameters  ------------------------------*/


int SVMgetB(double* b);
int SVMgetSvCount(int* svCount);
int SVMgetFeatureCount(int* featureCount);

/* In both functions below, the array is assumed to be properly allocated.
   Use SVMgetSvCount to get the correct length for alpha and svIndex arrays */

int SVMgetAlphas(SVMreal* alpha);
int SVMgetSvIndex(int* svIndex);



/*------------------------------------  Input/output  -------------------------------------*/

int SVMfreadClassifier(FILE *f, int binary);
int SVMfwriteClassifier(FILE *f, int binary);

int SVMreadClassifier(char* fileName, int binary);
int SVMwriteClassifier(char* fileName, int binary);


/*------------------------------------  Command line options  -------------------------------*/


void SVMprintSvmOptions();
void SVMparseSvmOptions(char* argv[], int argc);

int SVMgetBinaryFlag();

void SVMprintDataOptions();
void SVMprintDataOptionHelp();
int SVMparseDataOptions(char *argv[], int argc, int* k, int posCount, int negCount);

int SVMreadData(SVMreal** data, int rows, int cols);


#endif // __SVM_LIB_C_H__



