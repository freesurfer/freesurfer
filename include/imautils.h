#ifndef IMAUTILS_H
#define IMAUTILS_H


#define IMA_TYPE_SHORT  0
#define IMA_TYPE_INT    1
#define IMA_TYPE_LONG   2
#define IMA_TYPE_FLOAT  3
#define IMA_TYPE_DOUBLE 4
#define IMA_TYPE_STRING 5

typedef struct {
  char *FileName;
  char *PatientName;
  char *PatientDOB;
  char *StudyDate;
  char *StudyTime;
  char *PulseSequence;

  float FlipAngle;
  float EchoTime;
  float RepetitionTime;
  float InversionTime;

  int  StudyNo;
  int  SeriesNo;
  int  NFilesInSeries;
  int  ImageNo;    /* within the study, not series */
  int  NImageRows;
  int  NImageCols;
  float ImgCenter[3];
  float SliceThickness;
  float DistanceFactor;

  /* Volume/Run Related parameters */
  float Vc[3];     /* RAS col direction cosines*/
  float Vr[3];     /* RAS row direction cosines*/
  float Vs[3];     /* RAS slice direction cosines*/

  int   IsMosaic;    /* Image is a mosaic of slices */
  int   VolDim[3];   /* number of cols rows slices */
  float VolRes[3];   /* Resolution of col, row, slice in mm */
  int   NFrames;     /* Equals lRepetitions + 1 */

  //float VolCenter[3]; /* Exact RAS center of the volume */
  
  int   ErrorFlag;   /* Set for error, eg, aborted run */

} IMAFILEINFO;

/*******************************************/
typedef struct {
  char *key;
  int   offset;
  char *typestring;
  int   type;
  int   typesize;
  int   nitems;
} IMA_DICTIONARY_ENTRY;
#define NMAX_IMA_DICTIONARY 500

#ifdef IMAUTILS_SRC
  IMA_DICTIONARY_ENTRY ImaDictionary[NMAX_IMA_DICTIONARY];
  int nImaDictionary = 0, ImaDictionaryGood = 0;
  char *imaTypeString[6] = {"short","int","long","float","double","string"};
  int imaTypeSize[6] = {sizeof(short),sizeof(int),sizeof(long),
              sizeof(float),sizeof(double),sizeof(char)};
#else
  extern IMA_DICTIONARY_ENTRY ImaDictionary[NMAX_IMA_DICTIONARY];
  extern int   nImaDictionary, ImaDictionaryGood;
  extern char *imaTypeString[6];
  extern int   imaTypeSize[6];
#endif


int imaTypeFromString(char *typestring);
void *imaLoadVal(FILE *imafp, int offset, int nbytes, int nitems, void *pval);

void MkImaDictionary(void);
void DumpImaDictionary(FILE *fp);
int  DumpImaDictionaryVal(FILE *fp, char *imafile);
void *imaLoadValFromKey(FILE *imafp, char *key, void *pval);
int imaPrintVal(FILE *fp, int type, void *pval);
int imaTypeFromKey(char *key);

int imaIsSiemensIMA(char *imafile);
int imaParseName(char *imafile, int *StudyNo, int *SeriesNo, int *ImageNo);
int imaHasIMAExtension(char *filename);
int imaCountFilesInSeries(char *imafile, int *FirstImageNo);

IMAFILEINFO *imaLoadFileInfo(char *imafile);
int imaDumpFileInfo(FILE *fp, IMAFILEINFO *ifi);

#if 0 /*------------------------------------------------*/

typedef struct {
  int offset;
  int type;
  int nbytes;
  int nitems;
  int valgood;
  short  sval;
  int    ival;
  long   lval;
  float  fval;
  double dval;
  char  *string;
  char  *descr;
} IMAELEMENT;

typedef struct {
  char *filename;
  int nelements;
  IMAELEMENT **e;
} IMAFILEINFO;

int imaDumpElement(FILE *fp, IMAELEMENT *e);
int imaPrintElementValue(FILE *fp, IMAELEMENT *e);
int imaFreeElement(IMAELEMENT **ppe);
IMAELEMENT *imaMakeElement(char *descr, int offset, int type, int nitems);
IMAELEMENT *imaCopyElement(IMAELEMENT *esrc);
IMAELEMENT *imaLoadElement(FILE *imafp, int offset, int type, char *descr,
         int nitems);
int imaLoadElementVal(FILE *imafp, IMAELEMENT *e);
int imaAddElement(IMAFILEINFO *ifi, IMAELEMENT *e);
int imaDumpFileInfo(FILE *fp, IMAFILEINFO *ifi);
IMAFILEINFO *imaDefaultFileInfo(void);
IMAFILEINFO *imaLoadDefaultFileInfo(char *imafile);
#endif

#endif /*#ifndef IMAUTILS_H*/
