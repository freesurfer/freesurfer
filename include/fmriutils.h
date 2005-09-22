
#ifndef FMRIUTILS_INC
#define FMRIUTILS_INC

#include "matrix.h"
#include "mri.h"
#include "fsglm.h"

#ifdef X
#undef X
#endif

/*---------------------------------------------------------*/
typedef struct{
  GLMMAT *glm;       // Holds all the glm stuff
  MRI *y;            // Input data
  MATRIX *Xg;        // Global regressor matrix
  int npvr;          // Number of per-voxel regressors
  MRI *pvr[50];      // Per-voxel regressors (local)
  float DOF;         // DOF
  int Xcols;         // X->cols + npvr
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

  int ncontrasts;    // Number of contrasts
  char *cname[100];  // Contrast names
  MATRIX *C[100];    // Contrast matrices
  MRI *gamma[100];   // gamma = C*beta
  MRI *F[100];       // F = gamma'*iXtX*gamma/(rvar*J)
  MRI *sig[100];     // sig = significance of the F
} MRIGLM;
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
int MRIglmFit(MRIGLM *glmmri);

#endif
