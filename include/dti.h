// $Id: dti.h,v 1.2 2006/09/08 23:10:42 greve Exp $

#ifndef DTI_INC
#define DTI_INC

typedef struct {
  float bValue;
  int nAcq;
  int nDir;
  int DiffMode;
  char *GradFile;
  MATRIX *GradDir;
  MATRIX *GradDirNorm;
  MATRIX *B;
} DTI;

const char *DTIsrcVersion(void);
int DTIparamsFromSiemensAscii(char *fname, float *bValue, 
			      int *nAcq, int *nDir, int *DiffMode);
int DTIloadGradients(DTI *dti, char *GradFile);
DTI *DTIstructFromSiemensAscii(char *fname);
int DTInormGradDir(DTI *dti);
int DTIdesignMatrix(DTI *dti);
#endif //#ifndef FSENV_INC
