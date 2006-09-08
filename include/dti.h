// $Id: dti.h,v 1.1 2006/09/08 22:21:26 greve Exp $

#ifndef DTI_INC
#define DTI_INC

const char *DTIsrcVersion(void);
int DTIparamsFromSiemensAscii(char *fname, float *bValue, 
			      int *nAcq, int *nDir, int *DiffMode);


#endif //#ifndef FSENV_INC
