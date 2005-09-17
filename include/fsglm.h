// fsglm.h - include file for fsglm.c
// $Id: fsglm.h,v 1.1 2005/09/17 23:10:27 greve Exp $

#ifndef FSGLM_H
#define FSGLM_H

#include "matrix.h"

const char * GLMSrcVersion(void);
#undef X

#define GLMMAT_NCONTRASTS_MAX 100
typedef struct{
  MATRIX *y;
  MATRIX *X;
  double DOF;
  int ill_cond_flag;

  MATRIX *beta;
  MATRIX *yhat;
  MATRIX *eres;
  double rvar;

  int ncontrasts;    // Number of contrasts
  MATRIX *C[GLMMAT_NCONTRASTS_MAX];    // Contrast matrices
  MATRIX *Ct[GLMMAT_NCONTRASTS_MAX];   // transposes Contrast matrices

  MATRIX *gamma[GLMMAT_NCONTRASTS_MAX];
  double F[GLMMAT_NCONTRASTS_MAX];
  double p[GLMMAT_NCONTRASTS_MAX];

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



