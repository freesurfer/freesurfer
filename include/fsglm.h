/**
 * @file  fsglm.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: greve $
 *    $Date: 2011/05/05 15:28:03 $
 *    $Revision: 1.15 $
 *
 * Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
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
// $Id: fsglm.h,v 1.15 2011/05/05 15:28:03 greve Exp $

#ifndef FSGLM_H
#define FSGLM_H

#include "matrix.h"

const char * GLMSrcVersion(void);
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
  char *Cname[GLMMAT_NCONTRASTS_MAX];    // Contrast names
  double Ccond[GLMMAT_NCONTRASTS_MAX];    // C condition number
  int UseGamma0[GLMMAT_NCONTRASTS_MAX];  // Flag
  MATRIX *gamma0[GLMMAT_NCONTRASTS_MAX];  // Expected value of gamma

  int ypmfflag[GLMMAT_NCONTRASTS_MAX];    // flag to compute PMF
  MATRIX *Mpmf[GLMMAT_NCONTRASTS_MAX];    // Contrast PMF matrices
  MATRIX *ypmf[GLMMAT_NCONTRASTS_MAX];

  MATRIX *gamma[GLMMAT_NCONTRASTS_MAX];
  double F[GLMMAT_NCONTRASTS_MAX];
  double p[GLMMAT_NCONTRASTS_MAX];

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
int GLMdump(char *dumpdir, GLMMAT *glm);
int GLMresynthTest(int niters, double *prvar);
MATRIX *GLMpmfMatrix(MATRIX *C, double *cond, MATRIX *P);
int GLMdof(GLMMAT *glm);



#endif



