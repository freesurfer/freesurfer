// fsglm.h - include file for fsglm.c
// $Id: fsglm.h,v 1.2 2005/09/18 04:54:35 greve Exp $

#ifndef FSGLM_H
#define FSGLM_H

#include "matrix.h"

const char * GLMSrcVersion(void);
#undef X

#define GLMMAT_NCONTRASTS_MAX 100
typedef struct{
  MATRIX *y;   // input: nframes-by-1 (can only be 1)
  MATRIX *X;   // Design matrix: nframes-by-ncols
  // Note: weighted GLM not included here. To do weighted,
  // weight y and X prior to analysis.
  double DOF;  // nrows of X - ncols
  int ill_cond_flag;

  MATRIX *beta;
  MATRIX *yhat;
  MATRIX *eres;
  double rvar;  // Residual error variance

  int ncontrasts;    // Number of contrasts
  MATRIX *C[GLMMAT_NCONTRASTS_MAX];    // Contrast matrices

  MATRIX *gamma[GLMMAT_NCONTRASTS_MAX];
  double F[GLMMAT_NCONTRASTS_MAX];
  double p[GLMMAT_NCONTRASTS_MAX];

  // These are matrices to hold intermediate values
  MATRIX *Ct[GLMMAT_NCONTRASTS_MAX];   // transposes of contrast matrices
  MATRIX *Xt,*XtX,*iXtX,*Xty;
  MATRIX *CiXtX[GLMMAT_NCONTRASTS_MAX];
  MATRIX *CiXtXCt[GLMMAT_NCONTRASTS_MAX];

  MATRIX *gammat[GLMMAT_NCONTRASTS_MAX];
  MATRIX *gCVM[GLMMAT_NCONTRASTS_MAX];
  MATRIX *igCVM[GLMMAT_NCONTRASTS_MAX];
  MATRIX *gtigCVM[GLMMAT_NCONTRASTS_MAX];

} GLMMAT;

GLMMAT *GLMalloc(void);
int GLMfree(GLMMAT **pgm);
int GLMtransposeC(GLMMAT *gm);
int GLMfit(GLMMAT *glm);
int GLMtest(GLMMAT *glm);
int GLMprofile(int nrows, int ncols, int ncon, int niters);

GLMMAT *GLMsynth(void);
int GLMdump(char *dumpdir, GLMMAT *glm);
int GLMresynthTest(int niters, double *prvar);

#endif



