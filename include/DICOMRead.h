/**
 * @brief DICOM 3.0 reading functions
 *
 */
/*
 * Original Author: Sebastien Gicquel and Douglas Greve, 06/04/2001
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


#ifndef _DICOMRead_H
#define _DICOMRead_H


#include "dicom.h"
#include "lst.h"
#include "dicom_objects.h"
#include "condition.h"

#include "nii_dicom.h"
#include "nii_dicom_batch.h"

#define NUMBEROFTAGS 24
#define SHORTSIZE 16
#ifndef INTSIZE
#define INTSIZE 16
#endif
#define LONGSIZE 32

#define DCM_NOERROR 0
#define DCM_NODICOMFILES 1
#define DCM_MULTIPLESTUDIES 2

typedef unsigned short int BOOL;
// typedef unsigned short int bool;

// BEVIN
#ifndef true

// These two lines are the original lines
#define true 1
#define false 0

#elif true != 1 || false != 0

#error "Incompatible definitions of true or false"

#endif
// END OF BEVIN'S CHANGES

#ifdef _DICOMRead_SRC
char *SDCMStatusFile = 0;
char *SDCMListFile = 0;
int  UseDICOMRead2 = 1; // use new dicom reader by default
int  UseDCM2NIIX = 1; // changed to 1 on 4/06/2023
const char *DCM2NIIX_outdir = NULL;
int  DCM2NIIX_createBIDS = 0;
/* These variables allow the user to change the first tag checked to
   get the slice thickness.  This is needed with siemens mag res
   angiogram (MRAs) */
long SliceResElTag1 = 0x88; //Spacing Between Slices, only makes sense for 2D multislice or 3D multislab
long SliceResElTag2 = 0x50; //Slice Thickness, make sense for 3D
int AutoSliceResElTag = 0; // automatically determine which tag to use based on 18,23
#else
extern char *SDCMStatusFile;
extern char *SDCMListFile;
extern int  UseDICOMRead2;
extern int  UseDCM2NIIX;
extern const char *DCM2NIIX_outdir;
extern int  DCM2NIIX_createBIDS;
extern long SliceResElTag1;
extern long SliceResElTag2;
extern int AutoSliceResElTag;
#endif

typedef enum
{
  // note: increment the #DEFINE NUMBEROFTAGS at the top of this file
  // if adding elements to this DCM_TagList

  // general infos
  DCM_StudyDate,
  DCM_PatientName,
  DCM_Manufacturer,
  DCM_StudyTime,
  DCM_SeriesTime,
  DCM_AcquisitionTime,

  // image data format identifier (raw, or JPEG compressed)
  DCM_TransferSyntaxUID,

  // image dimensions
  DCM_SliceThickness,
  DCM_xsize,
  DCM_ysize,
  DCM_ImageNumber,
  DCM_SeriesNumber,
  DCM_Rows,
  DCM_Columns,
  DCM_BitsAllocated,
  DCM_NumberOfFrames,
  DCM_FieldOfView,

  // image position and orientation
  DCM_ImagePosition,
  DCM_ImageOrientation,

  // acquisition parameters
  DCM_EchoTime,
  DCM_RepetitionTime,
  DCM_InversionTime,
  DCM_FlipAngle,
  DCM_EchoNumber

}
DCM_TagList;

// This structure is for generic dicoms (see below for siemens specific)
typedef struct
{
  // DICOM file name
  char *FileName;

  // general infos
  char *StudyDate,
  *PatientName,
  *Manufacturer,
  *StudyTime,
  *SeriesTime,
  *AcquisitionTime;

  // image data format identifier (raw, or JPEG compressed)
  char *TransferSyntaxUID;

  // image dimensions
  double SliceThickness,
  xsize,
  ysize,
  FieldOfView;
  unsigned short ImageNumber,
  Rows,
  Columns,
  BitsAllocated,
  NumberOfFrames,
  SeriesNumber;

  // image position and orientation
  double ImagePosition[3],
  ImageOrientation[6],
  FirstImagePosition[3],
  LastImagePosition[3],
  Vc[3],Vr[3],Vs[3];

  // acquisition parameters
  double EchoTime,
    RepetitionTime,
    InversionTime,
    FlipAngle,
    FieldStrength;
  short EchoNumber;
  char *PhEncDir;
  double bval, bvecx, bvecy, bvecz;

  // pixels
  void *PixelData;
  unsigned char min8,  max8;
  unsigned short int min16, max16;

  // Rescaling parameters
  double RescaleIntercept, RescaleSlope; //(0028,1052) (0028,1053)

}
DICOMInfo ;

/*--- Relevant Info for a single Siemens DICOM File ------*/
typedef struct
{
  char *FileName;
  char *PatientName;
  char *StudyDate;
  char *StudyTime;
  char *SeriesTime;
  char *AcquisitionTime;
  char *PulseSequence;
  char *ProtocolName;
  char *PhEncDir;
  char *NumarisVer;
  char *ScannerModel;

  // This stores the 'Transfer Syntax Unique Identification',
  // which reports the structure of the image data, revealing
  // whether the data has been compressed. See:
  // http://www.psychology.nottingham.ac.uk/staff/cr1/dicom.html
  char *TransferSyntaxUID;

  int   EchoNo;
  float FlipAngle;
  float EchoTime;
  float RepetitionTime;
  float InversionTime;
  float FieldStrength;

  float PhEncFOV;
  float ReadoutFOV;

  int  SeriesNo;
  int  ImageNo;
  int  NImageRows;
  int  NImageCols;
  float ImgPos[3];

  int  lRepetitions;
  int  SliceArraylSize;

  /* Volume/Run Related parameters */
  float Vc[3];     /* RAS col direction cosines*/
  float Vr[3];     /* RAS row direction cosines*/
  float Vs[3];     /* RAS slice direction cosines*/

  int   RunNo;       /* Run Number that this is associated with */
  int   IsMosaic;    /* Image is a mosaic of slices */
  int   VolDim[3];   /* number of cols rows slices */
  float VolRes[3];   /* Resolution of col, row, slice in mm */
  float VolCenter[3]; /* Exact RAS center of the volume */
  int   NFrames;     /* Equals lRepetitions + 1 */
  double bValue;        /* for DWI */
  int    nthDirection;  /* also for DWI */
  int UseSliceScaleFactor; /* for slice-by-slice scaling (0020,4000) */
  double SliceScaleFactor; /* for slice-by-slice scaling (0020,4000) */
  double bval, bvecx, bvecy, bvecz;
  float LargestValue; // 0x28, 0x107
  int   ErrorFlag;   /* Set for error, eg, aborted run */

  // Rescaling parameters
  double RescaleIntercept, RescaleSlope; //(0028,1052) (0028,1053)

}
SDCMFILEINFO;

#ifdef _DICOMRead_SRC
char DICOMReadFirstDicomFile[5000];
#else
extern char DICOMReadFirstDicomFile[5000];
#endif

void PrintDICOMInfo(DICOMInfo *dcminfo);
CONDITION GetString(DCM_OBJECT** object, DCM_TAG tag, char **st);
CONDITION GetUSFromString(DCM_OBJECT** object, 
                          DCM_TAG tag, 
                          unsigned short *us);
CONDITION GetShortFromString(DCM_OBJECT** object, DCM_TAG tag, short *sh);
CONDITION GetUSFromUS(DCM_OBJECT** object, DCM_TAG tag, unsigned short *us);
CONDITION GetShortFromShort(DCM_OBJECT** object, DCM_TAG tag, short *ss);
CONDITION GetPixelData_Save(DCM_OBJECT** object, 
                            DCM_TAG tag, 
                            unsigned short **ad);
CONDITION GetPixelData(DCM_OBJECT** object, DCM_TAG tag, void **ad);
CONDITION GetDoubleFromString(DCM_OBJECT** object, DCM_TAG tag, double *d);
CONDITION GetMultiDoubleFromString(DCM_OBJECT** object, 
                                   DCM_TAG tag, 
                                   double *d[], 
                                   int multiplicity);
CONDITION GetMultiShortFromString(DCM_OBJECT** object, 
                                  DCM_TAG tag, 
                                  short *us[], 
                                  int multiplicity);
CONDITION GetDICOMInfo(const char *fname, 
                       DICOMInfo *dcminfo, 
                       BOOL ReadImage, 
                       int ImageNumber);
void *ReadDICOMImage(int nfiles, DICOMInfo **aDicomInfo);
void SortFiles(char *fNames[], 
               int nFiles,
               DICOMInfo ***ptrDicomArray, 
               int *nStudies);
int IsDICOM(const char *fname);
int ScanDir(const char *PathName, char ***FileNames, int *NumberOfFiles);
int CleanFileNames(char **FileNames, 
                   int NumberOfDICOMFiles, 
                   char ***CleanedFileNames);
int DICOMRead(const char *FileName, MRI **mri, int ReadImage);

int SortDCMFileInfo(DICOMInfo **dcmfi_list, int nlist);
int CompareDCMFileInfo(const void *a, const void *b);
int DCMCountFrames(DICOMInfo **dcmfi_list, int nlist);
int DCMSliceDir(DICOMInfo **dcmfi_list, int nlist);
MRI *DICOMRead2(const char *dcmfile, int LoadVolume);
MRIFSSTRUCT *DICOMRead3(const char *dcmfile, int LoadVolume);

DCM_ELEMENT *GetElementFromFile(const char *dicomfile, long grpid, long elid);
int AllocElementData(DCM_ELEMENT *e);
char *ElementValueString(DCM_ELEMENT *e, int DoBackslash);
int FreeElementData(DCM_ELEMENT *e);
DCM_ELEMENT *GetElementFromFile(const char *dicomfile, long grpid, long elid);
DCM_OBJECT *GetObjectFromFile(const char *fname, unsigned long options);
int IsSiemensDICOM(const char *dcmfile);
char *SiemensAsciiTag(const char *dcmfile,const  char *TagString, int flag);
char *SiemensAsciiTagEx(const char *dcmfile,const  char *TagString, int cleanup);
int dcmGetNCols(const char *dcmfile);
int dcmGetNRows(const char *dcmfile);
int dcmGetVolRes(const char *dcmfile, float *ColRes, float *RowRes, float *SliceRes);
int dcmImageDirCos(const char *dcmfile,
                   float *Vcx, float *Vcy, float *Vcz,
                   float *Vrx, float *Vry, float *Vrz);
int sdcmSliceDirCos(const char *dcmfile, float *Vsx, float *Vsy, float *Vsz);
int dcmImagePosition(const char *dcmfile, float *x, float *y, float *z);
int sdcmIsMosaic(const char *dcmfile, 
                 int *pNcols, 
                 int *pNrows, 
                 int *pNslices, 
                 int *pNframes);
MATRIX *sdcmAutoAlignMatrix(const char *dcmfile);

int DumpSDCMFileInfo(FILE *fp, SDCMFILEINFO *sdcmfi);
int FreeSDCMFileInfo(SDCMFILEINFO **ppsdcmfi);
SDCMFILEINFO *GetSDCMFileInfo(const char *dcmfile);
SDCMFILEINFO **ScanSiemensDCMDir(const char *PathName, int *NSDCMFiles);
int CompareSDCMFileInfo(const void *a, const void *b);
int SortSDCMFileInfo(SDCMFILEINFO **sdcmfi_list, int nlist);

int sdfiAssignRunNo(SDCMFILEINFO **sdcmfi_list, int nfiles);
int sdfiAssignRunNo2(SDCMFILEINFO **sdfi_list, int nlist);
int sdfiRunNo(const char *dcmfile, SDCMFILEINFO **sdfi_list, int nlist);
int sdfiNFilesInRun(const char *dcmfile, SDCMFILEINFO **sdfi_list, int nlist);
int sdfiCountFilesInRun(int RunNo, SDCMFILEINFO **sdfi_list, int nlist);
int *sdfiRunFileList(const char *dcmfile, SDCMFILEINFO **sdfi_list,
                     int nlist, int *NRunList);
MRI * sdcmLoadVolume(const char *dcmfile, int LoadVolume, int nthonly);
MRI * sdcmLoadVolumeAutoScale(const char *dcmfile, int LoadVolume, int nthonly);
int sdfiVolCenter(SDCMFILEINFO *sdfi);
int sdfiFixImagePosition(SDCMFILEINFO *sdfi);
int sdfiSameSlicePos(SDCMFILEINFO *sdfi1, SDCMFILEINFO *sdfi2);
int sdfiCountRuns(SDCMFILEINFO **sdfi_list, int nlist);
char *sdfiFirstFileInRun(int RunNo, SDCMFILEINFO **sdfi_list, int nlist);
int *sdfiRunNoList(SDCMFILEINFO **sdfi_list, int nlist, int *NRuns);
char **ScanSiemensSeries(const char *dcmfile, int *nList);
char **ReadSiemensSeries(const char *ListFile, int *nList,const char *dcmfile);
SDCMFILEINFO **LoadSiemensSeriesInfo(char **SeriesList, int nList);
char *sdcmExtractNumarisVer(const char *e_18_1020, int *Maj, int *Min, int *MinMin);
int sdfiIsSliceOrderReversed(SDCMFILEINFO *sdfi);
int dcmGetDWIParams(DCM_OBJECT *dcm, double *pbval, double *pxbvec, double *pybvec, double *pzbvec);
int dcmGetDWIParamsGE(DCM_OBJECT *dcm, double *pbval, double *pxbvec, double *pybvec, double *pzbvec);
int dcmGetDWIParamsPhilips(DCM_OBJECT *dcm, double *pbval, double *pxbvec, double *pybvec, double *pzbvec);
int dcmGetDWIParamsSiemens(DCM_OBJECT *dcm, double *pbval, double *pxbvec, double *pybvec, double *pzbvec);
int dcmGetDWIParamsSiemensAlt(DCM_OBJECT *dcm, double *pbval, double *pxbvec, double *pybvec, double *pzbvec);
int dcmImageDirCosObject(DCM_OBJECT *dcm, double *Vcx, double *Vcy, double *Vcz, double *Vrx, double *Vry, double *Vrz);
MATRIX *ImageDirCos2Slice(double Vcx, double Vcy, double Vcz,
			  double Vrx, double Vry, double Vrz,
			  double *Vsx, double *Vsy, double *Vsz);

int DCMcheckInterceptSlope(DCM_OBJECT *object);

#ifdef SunOS
/* kteich - this typedef is to keep the compiler from complaining
   about a struct defined in the parameter list... i don't know. i
   just replaced all 'struct dirent' in the code with 'struct_dirent'
   to reference the new typedef. seems to work. */
typedef struct dirent struct_dirent;

/* nicks - this code complains when built on Solaris 10.5
   it seems to already have scandir and alphasort (29 march 2006) */
#define HAVE_SCANDIR   1
#define HAVE_ALPHASORT 1

#ifndef HAVE_SCANDIR
int scandir(const char *dir, struct_dirent ***namelist,
            int (*select)(const struct_dirent *),
            int (*compar)(const void *, const void *));
#endif
#ifndef HAVE_ALPHASORT
int alphasort(const void *a, const void *b);
#endif

#endif

#endif
