/* fsgdf.h - header for freesurfer group descriptor file format */

#ifndef FSGDF_INC
#define FSGDF_INC

#ifdef X
#undef X
#endif

#include "matrix.h"

#define FSGDF_NCLASSES_MAX 100
#define FSGDF_NVARS_MAX    100
#define FSGDF_NINPUTS_MAX 1000

typedef struct {
  int version;
  char title[200];
  char measname[200];
  char regsubj[200];
  char datafile[1000];
  int  nclasses;
  char classlabel[FSGDF_NCLASSES_MAX][50]; /* [class][length]*/
  char classmarker[FSGDF_NCLASSES_MAX][50];  /* [class][length]*/  
  char classcolor[FSGDF_NCLASSES_MAX][50]; /* [class][length]*/  
  int  nvariables;
  char varlabel[FSGDF_NVARS_MAX][50]; /* [class][length]*/
  char defvarlabel[50]; /* default variable */
  int  ninputs;
  char subjid[FSGDF_NINPUTS_MAX][100];
  int  subjclassno[FSGDF_NINPUTS_MAX];
  float varvals[FSGDF_NINPUTS_MAX][FSGDF_NVARS_MAX];
} GROUPDESCRIPTOR, FSGD;

FSGD *gdfAlloc(int version);
int   gdfFree(FSGD **ppgd);
FSGD *gdfRead(char *gdfname);
int   gdfPrint(FILE *fp, FSGD *gd);
FSGD *gdfRead(char *gdfname);
MATRIX *gdfMatrixDOSS(FSGD *gd, MATRIX *X);
MATRIX *gdfMatrixDODS(FSGD *gd, MATRIX *X);

#endif //#ifndef FSGDF_INC


