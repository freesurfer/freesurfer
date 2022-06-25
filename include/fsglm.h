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


// fsglm.h - include file for fsglm.c

#ifndef FSGLM_H
#define FSGLM_H

#include "matrix.h"

#undef X

#define GLMMAT_NCONTRASTS_MAX 100
typedef struct {
  MATRIX *y;   // input: nframes-by-1 (can only be 1)
  MATRIX *X;   // Design matrix: nframes-by-ncols
  // Note: weighted GLM not included here. To do weighted,
  // weight y and X prior to analysis.
  double dof;  // nrows of X - ncols
  int AllowZeroDOF; 
  int ill_cond_flag;

  MATRIX *beta;
  MATRIX *yhat;
  MATRIX *eres;
  double rvar;  // Residual error variance

  MATRIX *yffxvar; // fixed-effects variance of each y
  int ffxdof;

  int ncontrasts;    // Number of contrasts
  MATRIX *C[GLMMAT_NCONTRASTS_MAX];    // Contrast matrices
  const char *Cname[GLMMAT_NCONTRASTS_MAX];    // Contrast names
  double Ccond[GLMMAT_NCONTRASTS_MAX];    // C condition number
  int UseGamma0[GLMMAT_NCONTRASTS_MAX];  // Flag
  MATRIX *gamma0[GLMMAT_NCONTRASTS_MAX];  // Expected value of gamma

  int ypmfflag[GLMMAT_NCONTRASTS_MAX];    // flag to compute PMF
  MATRIX *Mpmf[GLMMAT_NCONTRASTS_MAX];    // Contrast PMF matrices
  MATRIX *ypmf[GLMMAT_NCONTRASTS_MAX];

  MATRIX *gamma[GLMMAT_NCONTRASTS_MAX];
  double F[GLMMAT_NCONTRASTS_MAX];
  double p[GLMMAT_NCONTRASTS_MAX];
  double z[GLMMAT_NCONTRASTS_MAX]; // z derived from p
  double pcc[GLMMAT_NCONTRASTS_MAX]; // partial correlation coef

  /* When ReScaleX=1, rescale cols of X before computing inv(X'*X), 
     then rescale X'*X and inv(X'*X) so that it is transparent.*/
  int ReScaleX; 

  // These are matrices to hold intermediate values
  MATRIX *Ct[GLMMAT_NCONTRASTS_MAX];   // transposes of contrast matrices
  MATRIX *Xt,*XtX,*iXtX,*Xty;
  MATRIX *CiXtX[GLMMAT_NCONTRASTS_MAX];
  MATRIX *CiXtXCt[GLMMAT_NCONTRASTS_MAX];

  MATRIX *gammat[GLMMAT_NCONTRASTS_MAX];
  MATRIX *gCVM[GLMMAT_NCONTRASTS_MAX];
  MATRIX *igCVM[GLMMAT_NCONTRASTS_MAX];
  MATRIX *gtigCVM[GLMMAT_NCONTRASTS_MAX];

  // These are elements used to compute pcc
  int DoPCC;
  MATRIX *XCt[GLMMAT_NCONTRASTS_MAX]; // X*C'
  MATRIX *Dt[GLMMAT_NCONTRASTS_MAX]; // Null space of C
  MATRIX *XDt[GLMMAT_NCONTRASTS_MAX]; // X*D'
  MATRIX *RD[GLMMAT_NCONTRASTS_MAX]; // Residual forming matrix of XDt
  MATRIX *Xcd[GLMMAT_NCONTRASTS_MAX]; // RD * XCt
  MATRIX *Xcdt[GLMMAT_NCONTRASTS_MAX]; // Xcd'
  MATRIX *sumXcd[GLMMAT_NCONTRASTS_MAX]; // sum of columns of Xcd
  MATRIX *sumXcd2[GLMMAT_NCONTRASTS_MAX]; // sum of columns squared of Xcd
  MATRIX *yhatd[GLMMAT_NCONTRASTS_MAX]; // RD*yhat
  MATRIX *Xcdyhatd[GLMMAT_NCONTRASTS_MAX]; // Xcd*yhatd
  MATRIX *sumyhatd[GLMMAT_NCONTRASTS_MAX]; // sum(yhatd)
  MATRIX *sumyhatd2[GLMMAT_NCONTRASTS_MAX]; // sum(yhatd.^2)
  int debug=0;
}
GLMMAT;

GLMMAT *GLMalloc(void);
int GLMfree(GLMMAT **pgm);
int GLMallocX(GLMMAT *glm, int nrows, int ncols);
int GLMallocY(GLMMAT *glm);
int GLMallocYFFxVar(GLMMAT *glm);
int GLMcMatrices(GLMMAT *glm);
int GLMxMatrices(GLMMAT *glm);
int GLMfit(GLMMAT *glm);
int GLMtest(GLMMAT *glm);
int GLMtestFFx(GLMMAT *glm);

int GLManalyze(GLMMAT *glm);

int GLMprofile(int nrows, int ncols, int ncon, int niters);

GLMMAT *GLMsynth(void);
int GLMdump(const char *dumpdir, GLMMAT *glm);
int GLMresynthTest(int niters, double *prvar);
MATRIX *GLMpmfMatrix(MATRIX *C, double *cond, MATRIX *P);
int GLMdof(GLMMAT *glm);



#endif



