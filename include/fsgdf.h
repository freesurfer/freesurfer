/*
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


/* fsgdf.h - header for freesurfer group descriptor file format */

#ifndef FSGDF_INC
#define FSGDF_INC

#include "matrix.h"
#include "mri.h"

#ifdef X
#undef X
#endif

//This is for allowing repeats of subjects in the gdf for testing
//This will affect gdfCheckSubjRep() in fsgdf.c only.
#ifdef FSGDF_SRC
int fsgdf_AllowSubjRep = 0;
#else
extern int fsgdf_AllowSubjRep;
#endif


#define FSGDF_NCLASSES_MAX  128
#define FSGDF_NVARS_MAX     128
#define FSGDF_NINPUTS_MAX  20000

#define FSGD_FACTOR_DISCRETE 1
#define FSGD_FACTOR_CONTINUOUS 2

#define FSGD_FACTOR_CONTRAST_MAIN 1
#define FSGD_FACTOR_CONTRAST_INTERACTION 2

typedef struct 
{
  char name[200]; // Factor name
  int type; // discrete or continuous
  int nLevels; // number of levels for discrete factor
  char *Levels[100]; // level names for discrete factor
  int nFactors;// number of dfactor names for continous factor interactions
  char *FactorNames[100]; // dfactor names for continous factor interactions
}
FSGD_FACTOR;

typedef struct 
{
  char name[200];
  int type; // main, interaction
  int nFactors; // number of factors in the contrast
  char *FactorNames[100];
  FSGD_FACTOR *Factors[100];
}
FSGD_FACTOR_CONTRAST;


typedef struct
{
  int version;
  char title[200];
  char measname[200];
  char tessellation[20]; /* surface or volume */
  char regsubj[200];
  char datafile[1000];
  int  nclasses;
  char classlabel[FSGDF_NCLASSES_MAX][100]; /* [class][length]*/
  char classmarker[FSGDF_NCLASSES_MAX][100];  /* [class][length]*/
  char classcolor[FSGDF_NCLASSES_MAX][100]; /* [class][length]*/
  int  nvariables;
  char varlabel[FSGDF_NVARS_MAX][100]; /* [class][length]*/
  char defvarlabel[50]; /* default variable */

  int  nvarsfromfile;
  char *tablefile[100];
  char *varfield[100];
  int  fieldcol[100];
  int  datacol[100];

  int  ninputs;
  char subjid[FSGDF_NINPUTS_MAX][100];
  int  subjclassno[FSGDF_NINPUTS_MAX];
  float varvals[FSGDF_NINPUTS_MAX][FSGDF_NVARS_MAX];
  double NPerClass[FSGDF_NCLASSES_MAX];
  double VarMeans[FSGDF_NVARS_MAX];
  double VarStds[FSGDF_NVARS_MAX];
  double ClassVarMeans[FSGDF_NCLASSES_MAX][FSGDF_NVARS_MAX];
  char DesignMatFile[1000]; /* actual matlab4 mat file */
  char DesignMatMethod[100]; /* creation method */
  char gd2mtx_method[5];  //dods or doss
  MATRIX *X, *T; /* design matrix, T = inv(X'*X)*X' */
  MRI *data;
  double ResFWHM;
  int LogY; // indicates whether nat log of y was used
  int DeMean; // remove mean from continuous variables
  int ReScale; // divide continuous variables by stddev
  int nContrasts;
  char *ContrastName[50];
  int FContrastNSub[50];
  char **FContrastSub[50];
  int IsFContrast[50];
  MATRIX *C[50];

  FSGD_FACTOR *Factors[100];
  int nFactors;
  FSGD_FACTOR_CONTRAST *fc[100];
  int nFC;
  
}
GROUPDESCRIPTOR, FSGD;

FSGD   *gdfAlloc(int version);
int     gdfFree(FSGD **ppgd);
FSGD   *gdfRead(const char *gdfname, int LoadData);
int     gdfWrite(const char *gdfname, FSGD *gd);
MRI    *gdfReadDataInfo(const char *gdfname);
int     gdfPrintHeader(FILE *fp, FSGD *gd);
int     gdfCheckMatrixMethod(const char *gd2mtx_method);
int     gdfCheckNPerClass(FSGD *gd);
int     gdfPrint(FILE *fp, FSGD *gd);
int     gdfPrintStdout(FSGD *gd);
MATRIX *gdfMatrixDOSS(FSGD *gd, MATRIX *X);
MATRIX *gdfMatrixDODS(FSGD *gd, MATRIX *X);
MATRIX *gdfContrastDODS(FSGD *fsgd, float *wClass, float *wCovar);
MATRIX *gdfContrastDOSS(FSGD *fsgd, float *wClass, float *wCovar);
MATRIX *gdfMatrix(FSGD *gd, const char *gd2mtx_method, MATRIX *X);
int     gdfOffsetSlope(FSGD *gd, int classno, int varno,
                       int c, int r, int s, float *offset, float *slope);
int gdfCountItemsOnLine(FILE *fp);
int gdfCountItemsInString(const char *str);
char *gdfGetNthItemFromString(const char *str, const int nth);
int gdfClassNo(FSGD *gd, char *class_number);
int gdfGetVarLabelNo(FSGD *gd, char *LabelName);
int gdfStringIndex(char *str, char **list, int nlist);

int gdfGetTitle(FSGD *gd, char *title);
int gdfGetMeasurementName(FSGD *gd, char *name);
int gdfGetSubjectName(FSGD *gd, char *name);
double gdfGetFWHM(FSGD *gd);
int    gdfGetLogY(FSGD *gd);
int gdfGetDataFileName(FSGD *gd, char *filename);
int gdfGetNumClasses(FSGD *gd, int *nclasses);
int gdfGetNthClassLabel(FSGD *gd, int nclass, char *label);
int gdfGetNthClassMarker(FSGD *gd, int nclass, char *marker);
int gdfGetNthClassColor(FSGD *gd, int nclass, char *color);
int gdfGetNumVariables(FSGD *gd, int *nvariables);
int gdfGetNthVariableLabel(FSGD *gd, int nvariable, char *label);
int gdfGetDefaultVariable(FSGD *gd, char *label);
int gdfGetDefaultVariableIndex(FSGD *gd, int *nvariable);
int gdfGetNumSubjects(FSGD *gd, int *nsubjects);
int gdfGetNthSubjectID(FSGD *gd, int nsubject, char *id);
int gdfGetNthSubjectClass(FSGD *gd, int nsubject, int *class_number);
int gdfGetNthSubjectNthValue(FSGD *gd, int nsubject,
                             int nvariable, float *value);
int gdfGetNthSubjectMeasurement(FSGD *gd, int nsubject,
                                int x, int y, int z, float *value);

FSGD *gdfSubSet(FSGD *infsgd, int nClasses, char **ClassList,
                int nVars, char **VarList);
char **gdfCopySubjIdppc(FSGD *fsgd);
char *gdfGetSDataFromTable(char *tablefile, char *field,
                           int fieldcol, int datacol);
int gdfGetDDataFromTable(char *tablefile, char *field,
                         int fieldcol, int datacol, double *data);
int gdfVarMeans(FSGD *gd);
int gdfClassVarMeans(FSGD *gd);

#endif //#ifndef FSGDF_INC


