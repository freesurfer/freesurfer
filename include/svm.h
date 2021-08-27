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


#ifndef SVM_H
#define SVM_H

#include "const.h"

typedef struct
{
  int     ninputs ;        /* dimensionality of input space */
  int     ntraining ;      /* # of training vectors */
  char    class1_name[STRLEN] ;
  char    class2_name[STRLEN] ;
  double  *w ;            /* weights - ninputs long */
  double  threshold ;
  double  *alpha ;        /* Lagrange multipliers */
  double  C ;             /* 'weight' on allowing datapoints in margin */
  int     type ;          /* type of kernel used */
  int     nsupport ;      /* # of support vectors */
  float   *ysupport ;       /* outputs for support vectors */
  float   **xsupport ;      /* support vectors [nsupport][ninputs] */
  double  *asupport ;       /* alphas for support vectors */
  int     extra_args ;    /* stuff to write into file for user */
  char    **args ;
  double  sigma ;         /* only used if type == KERNEL_RBF */
}
SVM, SUPPORT_VECTOR_MACHINE ;


#define DEFAULT_SVM_C       1
#define DEFAULT_SVM_TOL     1e-7
#define DEFAULT_SVM_SIGMA   4

#define SVM_KERNEL_LINEAR       0
#define SVM_KERNEL_RBF          1
#define SVM_KERNEL_POLYNOMIAL   2


SVM   *SVMalloc(int ninputs, char *c1_name, char *c2_name) ;
int   SVMtrain(SVM *svm, float **x, float *y, int ntraining, double C,
               double tol, int max_iter);
int   SVMwrite(SVM *svm, char *fname) ;
SVM   *SVMread(char *fname) ;
int   SVMfree(SVM **psvm) ;
double SVMclassify(SVM *svm, float *x) ;
double SVMconstraint(SVM *svm, float *y, int ntraining) ;

#endif
