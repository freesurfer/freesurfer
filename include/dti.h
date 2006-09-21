// $Id: dti.h,v 1.3 2006/09/21 04:15:07 greve Exp $

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
#endif //#ifndef FSENV_INC
