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



#ifndef DTI_INC
#define DTI_INC

#include "mri.h"

typedef struct
{
  MATRIX *bValue;
  MATRIX *GradDir;
  MATRIX *GradDirNorm;
  MATRIX *B;
  // These strictly relate to Siemens/MGH DTI
  int nB0;
  int nDir;
  char *GradFile;
}
DTI;

int DTIfree(DTI **pdti);
int DTIparamsFromSiemensAscii(const char *fname, float *bValue,int *nDir, int *nB0);
int DTIloadGradients(DTI *dti,const  char *GradFile);
DTI *DTIstructFromSiemensAscii(const char *fname);
int DTInormGradDir(DTI *dti);
int DTIdesignMatrix(DTI *dti);
MRI *DTIbeta2Tensor(MRI *beta, MRI *mask, MRI *tensor);
MRI *DTIbeta2LowB(MRI *beta, MRI *mask, MRI *lowb);
MRI *DTIsynthDWI(MATRIX *X, MRI *beta, MRI *mask, MRI *synth);

int DTItensor2Eig(MRI *tensor, MRI *mask,   MRI **evals,
                  MRI **evec1, MRI **evec2, MRI **evec3);
MRI *DTIeigvals2FA(MRI *evals, MRI *mask, MRI *FA);
MRI *DTIeigvals2RA(MRI *evals, MRI *mask, MRI *RA);
MRI *DTIeigvals2VR(MRI *evals, MRI *mask, MRI *VR);
MRI *DTIradialDiffusivity(MRI *evals, MRI *mask, MRI *RD);

MRI *DTItensor2ADC(MRI *tensor, MRI *mask, MRI *adc);
int DTIsortEV(float *EigVals, MATRIX *EigVecs);
int DTIfslBValFile(DTI *dti,const  char *bvalfname);
int DTIfslBVecFile(DTI *dti,const  char *bvecfname);
MRI *DTIivc(MRI *evec, MRI *mask, MRI *ivc);
MATRIX *DTIloadBValues(const char *bvalfile);
MATRIX *DTIloadBVectors(const char *bvecfile);
int DTIwriteBVectors(MATRIX *bvecs,const char *bvecfile);
int DTIwriteBValues(MATRIX *bvals,const  char *bvalfile);
DTI *DTIstructFromBFiles(const char *bvalfile,const  char *bvecfile);
int DTIparsePulseSeqName(const char *pulseseq, double *bValue, int *nthDirection);
int DTIisFSLBVec(const char *fname);
int DTIbvecChangeSpace(MRI *vol, int desired_bvec_space);

#endif //#ifndef FSENV_INC
