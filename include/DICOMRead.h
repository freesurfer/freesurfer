#ifndef _DICOMRead_H
#define _DICOMRead_H

typedef enum {
  // general infos
  DCM_StudyDate, DCM_PatientName, DCM_Manufacturer, DCM_StudyTime, DCM_SeriesTime, DCM_AcquisitionTime,
  
  // image dimensions
  DCM_SliceThickness, DCM_xsize, DCM_ysize, DCM_ImageNumber, DCM_Rows, DCM_Columns, DCM_BitsAllocated, 
  DCM_NumberOfFrames, DCM_FieldOfView,

  // image position and orientation
  DCM_ImagePosition,  DCM_ImageOrientation,

  // acquisition parameters
  DCM_EchoTime, DCM_RepetitionTime, DCM_InversionTime, DCM_FlipAngle, DCM_EchoNumber
} DCM_TagList;

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
  
  // image dimensions
  double SliceThickness,
    xsize,
    ysize,
    FieldOfView;
  unsigned short ImageNumber,
    Rows, 
    Columns, 
    BitsAllocated, 
    NumberOfFrames;

  // image position and orientation
  double ImagePosition[3], 
    ImageOrientation[6],
    FirstImagePosition[3],
    LastImagePosition[3];

  // acquisition parameters
  double EchoTime, 
    RepetitionTime,
    InversionTime,
    FlipAngle;
  short EchoNumber;

  // pixels
  void *PixelData;
  unsigned char min8, 
    max8;
  unsigned short int min16, 
    max16;

}  DICOMInfo ;

#define NUMBEROFTAGS 22
#define SHORTSIZE 16
#define INTSIZE 16
#define LONGSIZE 32

#define DCM_NOERROR 0
#define DCM_NODICOMFILES 1
#define DCM_MULTIPLESTUDIES 2

typedef unsigned short int BOOL;
typedef unsigned short int bool;

#define true 1
#define false 0

bool IsTagPresent[NUMBEROFTAGS];

void PrintDICOMInfo(DICOMInfo *dcminfo);
CONDITION GetString(DCM_OBJECT** object, DCM_TAG tag, char **st);
CONDITION GetUSFromString(DCM_OBJECT** object, DCM_TAG tag, unsigned short *us);
CONDITION GetShortFromString(DCM_OBJECT** object, DCM_TAG tag, short *sh);
CONDITION GetUSFromUS(DCM_OBJECT** object, DCM_TAG tag, unsigned short *us);
CONDITION GetShortFromShort(DCM_OBJECT** object, DCM_TAG tag, short *ss);
CONDITION GetPixelData_Save(DCM_OBJECT** object, DCM_TAG tag, unsigned short **ad);
CONDITION GetPixelData(DCM_OBJECT** object, DCM_TAG tag, void **ad);
CONDITION GetDoubleFromString(DCM_OBJECT** object, DCM_TAG tag, double *d);
CONDITION GetMultiDoubleFromString(DCM_OBJECT** object, DCM_TAG tag, double *d[], int multiplicity);
CONDITION GetMultiShortFromString(DCM_OBJECT** object, DCM_TAG tag, short *us[], int multiplicity);
CONDITION GetDICOMInfo(char *fname, DICOMInfo *dcminfo, BOOL ReadImage, int ImageNumber);
void *ReadDICOMImage(int nfiles, DICOMInfo **aDicomInfo);
void SortFiles(char *fNames[], int nFiles, DICOMInfo ***ptrDicomArray, int *nStudies);
int IsDICOM(char *fname);
int ScanDir(char *PathName, char ***FileNames, int *NumberOfFiles);
int CleanFileNames(char **FileNames, int NumberOfDICOMFiles, char ***CleanedFileNames);
int DICOMRead(char *FileName, MRI **mri, int ReadImage);

#endif
