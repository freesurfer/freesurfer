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
  MATRIX *wg;        // Global weight vector
  int skipweight;    // Don't use weight even if w != NULL
  MRI *mask;         // Only proc within mask
  int n_ill_cond;    // Number of ill-conditioned voxels

  MRI *yffxvar;      // Fixed effects variance of each y
  int ffxdof;        // Fixed effects DOF

  MRI *cond;         // condition of X'*X
  int condsave;      // Flag to compute and save cond
  MRI *beta;         // beta = inv(X'*X)*X'*y
  MRI *yhat;         // yhat = X*beta
  int yhatsave;      // Flag to save yhat
  MRI *eres;         // eres = y - yhat
  MRI *rvar;         // rvar = sum(eres.^2)/DOF;

  MRI *gamma[100];   // gamma = C*beta
  MRI *gammaVar[100]; // gamma variance (t-tests only)
  MRI *F[100];       // F = gamma'*inv(C*inv(XtX)C')*gamma/(rvar*J)
  MRI *p[100];       // p = significance of the F
  MRI *z[100];       // z derived from p
  MRI *pcc[100];     // partial correlation coeff
  MRI *ypmf[100];    // partial model fit for each contrast
  MRI *FrameMask;    // Exclude a frame at a voxel if 0
}
MRIGLM;
/*---------------------------------------------------------*/

MRI *fMRImatrixMultiply(MRI *inmri, MATRIX *M, MRI *outmri);
MRI *fMRIcovariance(MRI *fmri, int Lag, float DOFAdjust, MRI *mask, MRI *covar);

MRI *fMRIcomputeT(MRI *ces, MATRIX *X, MATRIX *C, MRI *var, MRI *t);
MRI *fMRIcomputeF(MRI *ces, MATRIX *X, MATRIX *C, MRI *var, MRI *F);
MRI *fMRIsigT(MRI *t, float DOF, MRI *p);
MRI *fMRIsigF(MRI *F, float DOF1, float DOF2, MRI *sig);
MRI *fMRIsumSquare(MRI *fmri, int Update, MRI *sumsqr);
MRI *fMRInskip(MRI *inmri, int nskip, MRI *outmri);
MRI *fMRIndrop(MRI *inmri, int ndrop, MRI *outmri);
MRI *fMRIframe(MRI *inmri, int frame, MRI *outmri);
MRI *fMRIinsertFrame(MRI *srcmri, int srcframe, MRI *fmri, int frame);


MATRIX *MRItoMatrix(MRI *mri, int c, int r, int s,
                    int Mrows, int Mcols, MATRIX *M);
MATRIX *MRItoSymMatrix(MRI *mri, int c, int r, int s, MATRIX *M);
int MRIfromMatrix(MRI *mri, int c, int r, int s, MATRIX *M, MRI *FrameMask);
int MRIfromSymMatrix(MRI *mri, int c, int r, int s, MATRIX *M);
MRI *MRInormWeights(MRI *w, int sqrtFlag, int invFlag, MRI *mask, MRI *wn);

int MRIglmFitAndTest(MRIGLM *mriglm);
int MRIglmFit(MRIGLM *glmmri);
int MRIglmTest(MRIGLM *mriglm);
int MRIglmLoadVox(MRIGLM *mriglm, int c, int r, int s, int LoadBeta, GLMMAT *glm);
int MRIglmNRegTot(MRIGLM *mriglm);
VECTOR *MRItoVector(MRI *mri, int c, int r, int s, VECTOR *v);
int MRIsetSign(MRI *invol, MRI *signvol, int frame);
MRI *MRIvolMax(MRI *invol, MRI *out);
MRI *MRIvolMaxIndex(MRI *invol, int base, MRI *mask, MRI *out);
MRI *MRIvolMin(MRI *invol, MRI *out);
MRI *MRIconjunct(MRI *invol, MRI *out);
double MRIframeMax(MRI *vol, int frame, MRI *mask, int absflag,
                   int *cmax, int *rmax, int *smax);
MRI *MRIframeMean(MRI *vol, MRI *volmn);
MRI *MRIframeMedian(MRI *vol, MRI *volmn);
MRI *fMRIdetrend(MRI *y, MATRIX *X);
MRI *fMRItemporalAR1(MRI *fmri, float DOFAdjust, MRI *mask, MRI *ar1);
MRI *fMRIspatialAR1(MRI *src, MRI *mask, MRI *ar1);
MRI *fMRIspatialAR2(MRI *src, MRI *mask, MRI *ar2);
int fMRIspatialAR1Mean(MRI *ar1, MRI *mask, double *car1mn,
                       double *rar1mn,double *sar1mn);
int fMRIspatialAR2Mean(MRI *src, MRI *mask, double *car2mn,
                       double *rar2mn,double *sar2mn);
MRI *fMRIaddOffset(MRI *in, MRI *offset, MRI *mask, MRI *out);
MRI *fMRIsubSample(MRI *f, int Start, int Delta, int Stop, MRI *fsub);
MRI *fMRIexcludeFrames(MRI *f, int *ExcludeFrames, int nExclude, MRI *fex);
MRI *MRIframeSum(MRI *vol, MRI *volsum);
MRI *fMRItemporalGaussian(MRI *src, double gstdmsec, MRI *targ);
MRI *fMRIkurtosis(MRI *y, MRI *mask);
MRI *MRIpkurtosis(MRI *kvals, int dof, MRI *mask, int nsamples);

MATRIX *ASLinterpMatrix(int ntp);
MATRIX *fMRItoMatrix(MRI *fmri, MATRIX *M);
int fMRIfromMatrix(MATRIX *M, MRI *fmri);
MRI *fMRIspatialCorMatrix(MRI *fmri);
MRI *fMRIdistance(MRI *mri, MRI *mask);
MRI *fMRIcumSum(MRI *inmri, MRI *mask, MRI *outmri);
MRI *fMRIcumTrapZ(MRI *y, MATRIX *t, MRI *mask, MRI *yz);
MATRIX *HalfLife2Weight(double HalfLifeMin, MATRIX *tSec);
MRI *MRIframeNorm(MRI *src, MRI *mask, MRI *fnorm);
MRI *fMRIxcorr(MRI *v1, MRI *v2, MRI *mask, MRI *xcorr);
MRI *SpatialINorm(MRI *vol, MRI *mask, MRI *outvol);
MRI *fMRIspatialARN(MRI *src, MRI *mask, int N, MRI *arN);

#endif
