// QDEC Query-Design-Estimation-Contrast
// $Id: qdecutils.h,v 1.5 2005/11/09 00:53:31 greve Exp $

#ifndef QDEC_INC
#define QDEC_INC

#ifdef X
#undef X
#endif

#define QDEC_NLEVELS_MAX 10

#define QDEC_FACTOR_TYPE_CONTINUOUS 1
#define QDEC_FACTOR_TYPE_DISCRETE   2

/*---------------------------------------------------------*/
typedef struct{
  char *name;
  int type;
  int nLevels;
  char LevelNames[QDEC_NLEVELS_MAX][100];
  int nInputs;
  double  *cValues; // continuous value for each input
  char   **dValues; // or discrete value for each input
} QDEC_FACTOR, QDECF;
/*---------------------------------------------------------*/

/*---------------------------------------------------------*/
typedef struct{
  char *name;
  int nInputs;        // Number of subjects
  char **InputIdList; // Subject Ids
  int nFactors;
  QDEC_FACTOR *Factors;
} QDEC_DESIGN, QDECD;
/*---------------------------------------------------------*/

/*---------------------------------------------------------*/
typedef struct{
  char *name;         // Design name
  char *measure;      // thickness, sulc, curv
  char *hemi;         // lh or rh
  char *df1;          // name of discrete factor 1
  char *df2;          // name of discrete factor 2
  char *cf1;          // name of continuous factor 1
  char *cf2;          // name of continuous factor 2
  int  nsmooth;       // number of smooth steps
} QDEC_DESIGN_GUI, QDECDGUI;
/*---------------------------------------------------------*/




const char *QDECsrcVersion(void);
int QDECisContinuousFactor(QDECF *F);
int QDECisDiscreteFactor(QDECF *F);
int QDECallocFactor(QDECF *F);
int QDECallocDesign(QDECD *D);
QDECD *QDECsynthDesign(int nInputs, int nFactors);
const char *QDECfactorTypeName(QDECF *F);
int QDECdumpDesign(FILE *fp, QDECD *D);
int QDECnClasses(QDECD *D);
int QDECnContinuousFactors(QDECD *D);
int QDECnDiscreteFactors(QDECD *D);
int QDECnRegressors(QDECD *D);
int QDECfsgdFile(char *fsgdf, QDECD *D);
char *QDEClevels2ClassName(QDECD *D, int *nthlevels);
char *QDECinputClassName(QDECD *D, int nthInput);
//const char *QDECcheckDesign(QDECD *D);
QDECD *QDECloadTable(char *tablebase);
char **QDECdiscreteFactorNames(QDECD *D, int *ndf);
char **QDECcontinuousFactorNames(QDECD *D, int *ncf);
int QDECfactorNumber(QDECD *D, char *FactorName);
QDECF *QDECcloneFactor(QDECF *F, QDECF *Fsrc);
QDECD *QDECextractDesign(QDECD *DD, QDECDGUI *dgui);

#endif
