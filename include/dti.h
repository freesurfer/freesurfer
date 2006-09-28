// $Id: dti.h,v 1.4 2006/09/28 04:29:38 greve Exp $

#ifndef DTI_INC
#define DTI_INC

typedef struct {
  float bValue;
  int nAcq;
  int nB0;
  int nDir;
  char *GradFile;
  MATRIX *GradDir;
  MATRIX *GradDirNorm;
  MATRIX *B;
} DTI;

const char *DTIsrcVersion(void);
int DTIparamsFromSiemensAscii(char *fname, float *bValue, 
			      int *nAcq, int *nDir, int *nB0);
int DTIloadGradients(DTI *dti, char *GradFile);
DTI *DTIstructFromSiemensAscii(char *fname);
int DTInormGradDir(DTI *dti);
int DTIdesignMatrix(DTI *dti);
MRI *DTIbeta2Tensor(MRI *beta, MRI *mask, MRI *tensor);
MRI *DTIbeta2LowB(MRI *beta, MRI *mask, MRI *lowb);
int DTItensor2Eig(MRI *tensor, MRI *mask,   MRI **evals, 
		  MRI **evec1, MRI **evec2, MRI **evec3);
MRI *DTIeigvals2FA(MRI *evals, MRI *mask, MRI *FA);


#endif //#ifndef FSENV_INC
