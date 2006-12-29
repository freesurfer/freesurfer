/**
 * @file  imautils.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 02:08:59 $
 *    $Revision: 1.8 $
 *
 * Copyright (C) 2002-2007,
 * The General Hospital Corporation (Boston, MA). 
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 * Bug reports: analysis-bugs@nmr.mgh.harvard.edu
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
  char *key;
  int   offset;
  char *typestring;
  int   type;
  int   typesize;
  int   nitems;
}
IMA_DICTIONARY_ENTRY;
#define NMAX_IMA_DICTIONARY 500

extern IMA_DICTIONARY_ENTRY ImaDictionary[NMAX_IMA_DICTIONARY];
extern int   nImaDictionary, ImaDictionaryGood;
extern char *imaTypeString[6];
extern int   imaTypeSize[6];

int imaTypeFromString(char *typestring);
void *imaLoadVal(FILE *imafp, int offset, int nbytes, int nitems, void *pval);

void MkImaDictionary(void);
void DumpImaDictionary(FILE *fp);
int  DumpImaDictionaryVal(FILE *fp, char *imafile);
void *imaLoadValFromKey(FILE *imafp, char *key, void *pval);
int imaPrintVal(FILE *fp, int type, void *pval);
int imaTypeFromKey(char *key);

int imaIsSiemensIMA(char *imafile);
int imaParseName(char *imafile, int *StudyNo, int *SeriesNo, int *ImageNo,
                 char *Separator);
int imaHasIMAExtension(char *filename);
int imaCountFilesInSeries(char *imafile, int *FirstImageNo);

IMAFILEINFO *imaLoadFileInfo(char *imafile);
int imaDumpFileInfo(FILE *fp, IMAFILEINFO *ifi);
short *imaReadPixelData(IMAFILEINFO *ifi, short *PixelData);

#endif /*#ifndef IMAUTILS_H*/
