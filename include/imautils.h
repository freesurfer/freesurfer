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


#ifndef IMAUTILS_H
#define IMAUTILS_H


#define IMA_TYPE_SHORT  0
#define IMA_TYPE_INT    1
#define IMA_TYPE_LONG   2
#define IMA_TYPE_FLOAT  3
#define IMA_TYPE_DOUBLE 4
#define IMA_TYPE_STRING 5

typedef struct
{
  char *FileName;
  char *PatientName;
  char *PatientDOB;
  char  PatientGender[7];
  char *StudyDate;
  char *StudyTime;
  char *PulseSequence;

  float FlipAngle;
  float EchoTime;
  float RepetitionTime;
  float InversionTime;

  int  StudyNo;
  int  SeriesNo;
  int  NFilesInSeries;    /* Actual number */
  int  NFilesInSeriesExp; /* Expected */
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
  int   NFrames;

  int   NFilesPerFrame;

  //float VolCenter[3]; /* Exact RAS center of the volume */

  int   ErrorFlag;   /* Set for error, eg, aborted run */

}
IMAFILEINFO;

/*******************************************/
typedef struct
{
  const char *key;
  int   offset;
  const char *typestring;
  int   type;
  int   typesize;
  int   nitems;
}
IMA_DICTIONARY_ENTRY;
#define NMAX_IMA_DICTIONARY 500

extern IMA_DICTIONARY_ENTRY ImaDictionary[NMAX_IMA_DICTIONARY];
extern int   nImaDictionary, ImaDictionaryGood;
extern const char *imaTypeString[6];
extern int   imaTypeSize[6];

int imaTypeFromString(const char *typestring);
void *imaLoadVal(FILE *imafp, int offset, int nbytes, int nitems, void *pval);

void MkImaDictionary(void);
void DumpImaDictionary(FILE *fp);
int  DumpImaDictionaryVal(FILE *fp,const  char *imafile);
void *imaLoadValFromKey(FILE *imafp,const  char *key, void *pval);
int imaPrintVal(FILE *fp, int type, void *pval);
int imaTypeFromKey(const char *key);

int imaIsSiemensIMA(const char *imafile);
int imaParseName(const char *imafile, int *StudyNo, int *SeriesNo, int *ImageNo,
                 char *Separator);
int imaHasIMAExtension(const char *filename);
int imaCountFilesInSeries(const char *imafile, int *FirstImageNo);

IMAFILEINFO *imaLoadFileInfo(const char *imafile);
int imaDumpFileInfo(FILE *fp, IMAFILEINFO *ifi);
short *imaReadPixelData(IMAFILEINFO *ifi, short *PixelData);

#endif /*#ifndef IMAUTILS_H*/
