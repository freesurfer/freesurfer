/**
 * @file  fmriutils.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 02:08:59 $
 *    $Revision: 1.21 $
 *
 * Copyright (C) 2002-2007,
 * The General Hospital Corporation (Boston, MA). 
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 * Bug reports: analysis-bugs@nmr.mgh.harvard.edu
 *
 */



#ifndef FMRIUTILS_INC
#define FMRIUTILS_INC

#include "matrix.h"
#include "mri.h"
#include "fsglm.h"

#ifdef X
#undef X
#endif

/*---------------------------------------------------------*/
typedef struct
{
  GLMMAT *glm;       // Holds all the glm stuff
  MRI *y;            // Input data
  MATRIX *Xg;        // Global regressor matrix
  int npvr;          // Number of per-voxel regressors
  MRI *pvr[50];      // Per-voxel regressors (local)
  int nregtot;       // Total number of regressors

  int pervoxflag;    // 1 if X is per-voxel
  int XgLoaded;      // 1 if Xg has been loaded into glm->X

  MRI *w;            // Per-voxel, per-input weight
  int skipweight;    // Don't use weight even if w != NULL
  MRI *mask;         // Only proc within mask
  int n_ill_cond;    // Number of ill-conditioned voxels

  MRI *cond;         // condition of X'*X
  int condsave;      // Flag to compute and save cond
  MRI *beta;         // beta = inv(X'*X)*X'*y
  MRI *yhat;         // yhat = X*beta
  int yhatsave;      // Flag to save yhat
  MRI *eres;         // eres = y - yhat
  MRI *rvar;         // rvar = sum(eres.^2)/DOF;

  MRI *gamma[100];   // gamma = C*beta
  MRI *F[100];       // F = gamma'*inv(C*inv(XtX)C')*gamma/(rvar*J)
  MRI *p[100];       // p = significance of the F
  MRI *ypmf[100];    // partial model fit for each contrast
}
MRIGLM;
/*---------------------------------------------------------*/

const char *fMRISrcVersion(void);
MRI *fMRImatrixMultiply(MRI *inmri, MATRIX *M, MRI *outmri);
MRI *fMRIvariance(MRI *fmri, float DOF, int RmMean, MRI *var);
MRI *fMRIcomputeT(MRI *ces, MATRIX *X, MATRIX *C, MRI *var, MRI *t);
MRI *fMRIcomputeF(MRI *ces, MATRIX *X, MATRIX *C, MRI *var, MRI *F);
MRI *fMRIsigT(MRI *t, float DOF, MRI *p);
MRI *fMRIsigF(MRI *F, float DOF1, float DOF2, MRI *sig);
MRI *fMRIsumSquare(MRI *fmri, int Update, MRI *sumsqr);
MRI *fMRInskip(MRI *inmri, int nskip, MRI *outmri);
MRI *fMRIndrop(MRI *inmri, int ndrop, MRI *outmri);
MATRIX *MRItoMatrix(MRI *mri, int c, int r, int s,
                    int Mrows, int Mcols, MATRIX *M);
MATRIX *MRItoSymMatrix(MRI *mri, int c, int r, int s, MATRIX *M);
int MRIfromMatrix(MRI *mri, int c, int r, int s, MATRIX *M);
int MRIfromSymMatrix(MRI *mri, int c, int r, int s, MATRIX *M);
MRI *MRInormWeights(MRI *w, int sqrtFlag, int invFlag, MRI *mask, MRI *wn);

int MRIglmFitAndTest(MRIGLM *mriglm);
int MRIglmFit(MRIGLM *glmmri);
int MRIglmTest(MRIGLM *mriglm);
int MRIglmLoadVox(MRIGLM *mriglm, int c, int r, int s, int LoadBeta);
int MRIglmNRegTot(MRIGLM *mriglm);
VECTOR *MRItoVector(MRI *mri, int c, int r, int s, VECTOR *v);
int MRIsetSign(MRI *invol, MRI *signvol, int frame);
MRI *MRIvolMax(MRI *invol, MRI *out);
double MRIframeMax(MRI *vol, int frame, MRI *mask, int absflag,
                   int *cmax, int *rmax, int *smax);
MRI *MRIframeMean(MRI *vol, MRI *volmn);
MRI *fMRIdetrend(MRI *y, MATRIX *X);
MRI *fMRIspatialAR1(MRI *src, MRI *mask, MRI *ar1);
int fMRIspatialAR1Mean(MRI *src, MRI *mask, double *car1mn,
                       double *rar1mn,double *sar1mn);
MRI *fMRIaddOffset(MRI *in, MRI *offset, MRI *mask, MRI *out);

#endif
