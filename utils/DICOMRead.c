/*******************************************************
   DICOM 3.0 reading functions
   Author: Sebastien Gicquel
   Date: 06/04/2001
*******************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <sys/file.h>
#include <sys/types.h>
#include <time.h>
#include <sys/timeb.h>
#include <sys/time.h>
#include <ctype.h>
#include <dirent.h>
#include <malloc.h>
#include "dicom.h"
#include "lst.h"
#include "dicom_objects.h"
#include "mri.h"
#include "mri_identify.h"
#include "DICOMRead.h"

//#define _DEBUG



/*******************************************************
   PrintDICOMInfo
   Author: Sebastien Gicquel
   Date: 06/04/2001
   input: structure DICOMInfo
   output: prints specific DICOM fields
*******************************************************/


void PrintDICOMInfo(DICOMInfo *dcminfo)
{
  int i;
  char str[256];

  printf("-------------------------------------------------\n");
  printf("DICOM meta-header\n\n");

  printf("file name\t\t\t%s\n", dcminfo->FileName);

  printf("Date and time\n");
  if (IsTagPresent[DCM_StudyDate])
    sprintf(str, "%s", dcminfo->StudyDate);
  else
    strcpy(str, "not found");
  printf("\tstudy date\t\t%s\n", str);

  if (IsTagPresent[DCM_StudyTime])
    sprintf(str, "%s", dcminfo->StudyTime);
  else
    strcpy(str, "not found");
  printf("\tstudy time\t\t%s\n", str);

  if (IsTagPresent[DCM_SeriesTime])
    sprintf(str, "%s", dcminfo->SeriesTime);
  else
    strcpy(str, "not found");
  printf("\tseries time\t\t%s\n", str);

  if (IsTagPresent[DCM_AcquisitionTime])
    sprintf(str, "%s", dcminfo->AcquisitionTime);
  else
    strcpy(str, "not found");
  printf("\tacquisition time\t%s\n", str);

  printf("Identification\n");
  if (IsTagPresent[DCM_PatientName])
    sprintf(str, "%s", dcminfo->PatientName);
  else
    strcpy(str, "not found");
  printf("\tpatient name\t\t%s\n", str);

  if (IsTagPresent[DCM_Manufacturer])
    sprintf(str, "%s", dcminfo->Manufacturer);
  else
    strcpy(str, "not found");
  printf("\tmanufacturer\t\t%s\n", str);

  printf("Dimensions\n");
  if (IsTagPresent[DCM_Rows])
    sprintf(str, "%d", dcminfo->Rows);
  else
    strcpy(str, "not found");
  printf("\tnumber of rows\t\t%s\n", str);

  if (IsTagPresent[DCM_Columns])
    sprintf(str, "%d", dcminfo->Columns);
  else
    strcpy(str, "not found");
  printf("\tnumber of columns\t%s\n", str);

  printf("\tnumber of frames\t%d\n", dcminfo->NumberOfFrames);

  if (IsTagPresent[DCM_xsize])
    sprintf(str, "%g", dcminfo->xsize);
  else
    strcpy(str, "not found");
  printf("\tpixel width\t\t%s\n",  str);

  if (IsTagPresent[DCM_ysize])
    sprintf(str, "%g", dcminfo->ysize);
  else
    strcpy(str, "not found");
  printf("\tpixel height\t\t%s\n",  str);

  if (IsTagPresent[DCM_SliceThickness])
    sprintf(str, "%g", dcminfo->SliceThickness);
  else
    strcpy(str, "not found");
  printf("\tslice thickness\t\t%s\n",  str);

  printf("\tfield of view\t\t%g\n", dcminfo->FieldOfView);

  if (IsTagPresent[DCM_ImageNumber])
    sprintf(str, "%d", dcminfo->ImageNumber);
  else
    strcpy(str, "not found");
  printf("\timage number\t\t%s (might be not reliable)\n", str);

  printf("Acquisition parameters\n");
  if (IsTagPresent[DCM_EchoTime])
    sprintf(str, "%g", dcminfo->EchoTime);
  else
    strcpy(str, "not found");
  printf("\techo time\t\t%s\n",  str);

  if (IsTagPresent[DCM_RepetitionTime])
    sprintf(str, "%g", dcminfo->RepetitionTime);
  else
    strcpy(str, "not found");
  printf("\trepetition time\t\t%s\n", str);

  if (IsTagPresent[DCM_InversionTime])
    sprintf(str, "%g", dcminfo->InversionTime);
  else
    strcpy(str, "not found");
  printf("\tinversion time\t\t%s\n", str);

  if (IsTagPresent[DCM_EchoNumber])
    sprintf(str, "%d", dcminfo->EchoNumber);
  else
    strcpy(str, "not found");
  printf("\techo number\t\t%s\n", str);

  if (IsTagPresent[DCM_FlipAngle])
    sprintf(str, "%g", dcminfo->FlipAngle);
  else
    strcpy(str, "not found");
  printf("\tflip angle\t\t%s\n", str);

  if (IsTagPresent[DCM_BitsAllocated])
    sprintf(str, "%d", dcminfo->BitsAllocated);
  else
    strcpy(str, "not found");
  printf("\tbits allocated\t\t%s\n", str);

  printf("Spatial informations\n");
  
  printf("\tfirst image position\t");
  if (IsTagPresent[DCM_ImagePosition])
    {
      for (i=0; i<3; i++)
  {
    printf("%g ", dcminfo->FirstImagePosition[i]);
  }
      printf("\n");
    }
  else
    printf("notfound\n");

  printf("\tlast image position\t");
  if (IsTagPresent[DCM_ImagePosition])
    {
      for (i=0; i<3; i++)
  {
    printf("%g ", dcminfo->LastImagePosition[i]);
  }
      printf("\n");
    }
  else
    printf("notfound\n");

  printf("\timage orientation\t");
  if (IsTagPresent[DCM_ImageOrientation])
    {
      for (i=0; i<6; i++)
  {
    printf("%g ", dcminfo->ImageOrientation[i]);
  }
      printf("\n");
    }
  else
    printf("notfound\n");

  printf("-------------------------------------------------\n\n");

}


CONDITION GetString(DCM_OBJECT** object, DCM_TAG tag, char **st)
{
  DCM_ELEMENT attribute;
  CONDITION cond;
  void* ctx;

  attribute.tag=tag;
  cond=DCM_GetElement(object, tag, &attribute);
  if (cond != DCM_NORMAL)
    return cond;
  *st=(char*)calloc(attribute.length+1, sizeof(char));
  attribute.d.string=*st;
  ctx=NULL;
  cond=DCM_GetElementValue(object, &attribute, &attribute.length, &ctx);
  return cond;

}

CONDITION GetUSFromString(DCM_OBJECT** object, DCM_TAG tag, unsigned short *us)
{
  DCM_ELEMENT attribute;
  CONDITION cond;
  char* s;
  void* ctx;

  attribute.tag=tag;
  cond=DCM_GetElement(object, tag, &attribute);
  if (cond != DCM_NORMAL)
    return cond;
  s=(char*)calloc(attribute.length+1, sizeof(char));
  attribute.d.string=s;
  ctx=NULL;
  cond=DCM_GetElementValue(object, &attribute, &attribute.length, &ctx);
  *us=(unsigned short)atoi(s);
  free(s);
  return cond;
}

CONDITION GetShortFromString(DCM_OBJECT** object, DCM_TAG tag, short *sh)
{
  DCM_ELEMENT attribute;
  CONDITION cond;
  char* s;
  void* ctx;

  attribute.tag=tag;
  cond=DCM_GetElement(object, tag, &attribute);
  if (cond != DCM_NORMAL)
    return cond;
  s=(char*)calloc(attribute.length+1, sizeof(char));
  attribute.d.string=s;
  ctx=NULL;
  cond=DCM_GetElementValue(object, &attribute, &attribute.length, &ctx);
  *sh=(short)atoi(s);
  free(s);
  return cond;
}

CONDITION GetUSFromUS(DCM_OBJECT** object, DCM_TAG tag, unsigned short *us)
{
  DCM_ELEMENT attribute;
  CONDITION cond;
  void* ctx, *ot;

  attribute.tag=tag;
  cond=DCM_GetElement(object, tag, &attribute);
  if (cond != DCM_NORMAL)
    return cond;
  ot=(void*)malloc(attribute.length+1);
  attribute.d.ot=ot;
  ctx=NULL;
  cond=DCM_GetElementValue(object, &attribute, &attribute.length, &ctx);
  *us=*attribute.d.us;
  free(ot);
  return cond;
}

CONDITION GetShortFromShort(DCM_OBJECT** object, DCM_TAG tag, short *ss)
{
  DCM_ELEMENT attribute;
  CONDITION cond;
  void* ctx, *ot;

  attribute.tag=tag;
  cond=DCM_GetElement(object, tag, &attribute);
  if (cond != DCM_NORMAL)
    return cond;
  ot=(void*)malloc(attribute.length+1);
  attribute.d.ot=ot;
  ctx=NULL;
  cond=DCM_GetElementValue(object, &attribute, &attribute.length, &ctx);
  *ss=*attribute.d.ss;
  free(ot);
  return cond;
}

CONDITION GetPixelData_Save(DCM_OBJECT** object, DCM_TAG tag, unsigned short **ad)
{
  DCM_ELEMENT attribute;
  CONDITION cond;
  void *ctx, *ot;

  attribute.tag=tag;
  cond=DCM_GetElement(object, tag, &attribute);
  if (cond != DCM_NORMAL)
    return cond;
  ot=(void*)malloc(attribute.length+1);
  attribute.d.ot=ot;
  ctx=NULL;
  cond=DCM_GetElementValue(object, &attribute, &attribute.length, &ctx);
  *ad=(unsigned short *)ot;
  return cond;
}

CONDITION GetPixelData(DCM_OBJECT** object, DCM_TAG tag, void **ad)
{
  DCM_ELEMENT attribute;
  CONDITION cond;
  void *ctx, *ot;

  attribute.tag=tag;
  cond=DCM_GetElement(object, tag, &attribute);
  if (cond != DCM_NORMAL)
    return cond;
  ot=(void*)malloc(attribute.length+1);
  attribute.d.ot=ot;
  ctx=NULL;
  cond=DCM_GetElementValue(object, &attribute, &attribute.length, &ctx);
  *ad=ot;
  return cond;
}

CONDITION GetDoubleFromString(DCM_OBJECT** object, DCM_TAG tag, double *d)
{
  DCM_ELEMENT attribute;
  CONDITION cond;
  char* s;
  void* ctx;

  attribute.tag=tag;
  cond=DCM_GetElement(object, tag, &attribute);
  if (cond != DCM_NORMAL)
    return cond;
  s=(char*)calloc(attribute.length+1, sizeof(char));
  attribute.d.string=s;
  ctx=NULL;
  cond=DCM_GetElementValue(object, &attribute, &attribute.length, &ctx);
  *d=atof(s);
  free(s);
  return cond;
}

CONDITION GetMultiDoubleFromString(DCM_OBJECT** object, DCM_TAG tag, double *d[], int multiplicity)
{
  DCM_ELEMENT attribute;
  CONDITION cond;
  char *s, *ss;
  void* ctx;
  int i, j, mult;

  attribute.tag=tag;
  cond=DCM_GetElement(object, tag, &attribute);
  if (cond != DCM_NORMAL)
    return cond;
  s=(char*)calloc(attribute.length+1, sizeof(char));
  ss=(char*)calloc(attribute.length+1, sizeof(char));
  attribute.d.string=s;
  ctx=NULL;
  cond=DCM_GetElementValue(object, &attribute, &attribute.length, &ctx);


  i=0;
  j=0;
  mult=0;
  while (mult<multiplicity && i<attribute.length) {
    j=0;
    while (s[i]!='\\' && i<attribute.length)
      ss[j++]=s[i++];
    i++;
    ss[j]='\0';
    (*d)[mult++]=atof(ss);
  }
  free(s);
  free(ss);
  return cond;
}

CONDITION GetMultiShortFromString(DCM_OBJECT** object, DCM_TAG tag, short *us[], int multiplicity)
{
  DCM_ELEMENT attribute;
  CONDITION cond;
  char *s, *ss;
  void* ctx;
  int i, j, mult;

  attribute.tag=tag;
  cond=DCM_GetElement(object, tag, &attribute);
  if (cond != DCM_NORMAL)
    return cond;
  s=(char*)calloc(attribute.length+1, sizeof(char));
  ss=(char*)calloc(attribute.length+1, sizeof(char));
  attribute.d.string=s;
  ctx=NULL;
  cond=DCM_GetElementValue(object, &attribute, &attribute.length, &ctx);

  i=0;
  j=0;
  mult=0;
  while (mult<multiplicity && i<attribute.length) {
    j=0;
    while (s[i]!='\\' && i<attribute.length)
      ss[j++]=s[i++];
    i++;
    ss[j]='\0';
    (*us)[mult++]=atoi(ss);
  }
  free(s);
  free(ss);
  return cond;
}

/*******************************************************
   GetDICOMInfo
   Author: Sebastien Gicquel
   Date: 06/04/2001
   input: file name, structure DICOMInfo, boolean, image number
   output: fills in structure DICOMInfo with DICOM meta-header fields
*******************************************************/

CONDITION GetDICOMInfo(char *fname, DICOMInfo *dcminfo, BOOL ReadImage, int ImageNumber)
{
  DCM_OBJECT** object=(DCM_OBJECT**)calloc(1, sizeof(DCM_OBJECT*));
  DCM_TAG tag;
  CONDITION cond, cond2=DCM_NORMAL;
  double *tmp=(double*)calloc(10, sizeof(double));
  short *itmp=(short*)calloc(3, sizeof(short));

  int i;

  cond=DCM_OpenFile(fname, DCM_PART10FILE, object);
  if (cond != DCM_NORMAL)
    cond=DCM_OpenFile(fname, DCM_ORDERLITTLEENDIAN, object);
  if (cond != DCM_NORMAL)
    cond=DCM_OpenFile(fname, DCM_ORDERBIGENDIAN, object);
  if (cond != DCM_NORMAL)
    cond=DCM_OpenFile(fname, DCM_FORMATCONVERSION, object);
  if (cond != DCM_NORMAL)
    {
      printf("not DICOM images, sorry...\n");
      exit(1);
    }

  dcminfo->FileName=(char *)calloc(strlen(fname)+1, sizeof(char));
  strcpy(dcminfo->FileName, fname);
  dcminfo->FileName[strlen(fname)]='\0';

  // manufacturer
  tag=DCM_MAKETAG(0x8, 0x70);
  cond=GetString(object, tag, &dcminfo->Manufacturer);
  if (cond != DCM_NORMAL) {
    dcminfo->Manufacturer=(char *)calloc(256, sizeof(char));
    strcpy(dcminfo->Manufacturer, "UNKNOWN");
    cond2=cond;
#ifdef _DEBUG
    printf("WARNING: tag Manufacturer not found\n"); 
#endif  
  }
  else
    IsTagPresent[DCM_Manufacturer]=true;


  // patient name
  tag=DCM_MAKETAG(0x10, 0x10);
  cond=GetString(object, tag, &dcminfo->PatientName);
  if (cond != DCM_NORMAL) {
    dcminfo->PatientName=(char *)calloc(256, sizeof(char));
    strcpy(dcminfo->PatientName, "UNKNOWN");
    cond2=cond;
#ifdef _DEBUG
    printf("WARNING: tag PatientName not found\n"); 
#endif  
  }
  else
    IsTagPresent[DCM_PatientName]=true;

  // study date
  tag=DCM_MAKETAG(0x8, 0x20);
  cond=GetString(object, tag, &dcminfo->StudyDate);
  if (cond != DCM_NORMAL || strcmp(dcminfo->StudyDate, "00000000")==0) {
    dcminfo->StudyDate=(char *)calloc(256, sizeof(char));
    strcpy(dcminfo->StudyDate, "UNKNOWN");
    cond2=cond;
#ifdef _DEBUG
    printf("WARNING: tag StudyDate not found\n"); 
#endif  
  }
  else
    IsTagPresent[DCM_StudyDate]=true;

  // study time
  tag=DCM_MAKETAG(0x8, 0x30);
  cond=GetString(object, tag, &dcminfo->StudyTime);
  if (cond != DCM_NORMAL || strcmp(dcminfo->StudyTime, "00000000")==0) {
    dcminfo->StudyTime=(char *)calloc(256, sizeof(char));
    strcpy(dcminfo->StudyTime, "UNKNOWN");
    cond2=cond;
#ifdef _DEBUG
    printf("WARNING: tag StudyTime not found\n"); 
#endif  
  }
  else
    IsTagPresent[DCM_StudyTime]=true;

  // series time
  tag=DCM_MAKETAG(0x8, 0x31);
  cond=GetString(object, tag, &dcminfo->SeriesTime);
  if (cond != DCM_NORMAL || strcmp(dcminfo->SeriesTime, "00000000")==0) {
    dcminfo->SeriesTime=(char *)calloc(256, sizeof(char));
    strcpy(dcminfo->SeriesTime, "UNKNOWN");
    cond2=cond;
#ifdef _DEBUG
    printf("WARNING: tag SeriesTime not found\n"); 
#endif  
  }
  else
    IsTagPresent[DCM_SeriesTime]=true;

  // acquisition time
  tag=DCM_MAKETAG(0x8, 0x31);
  cond=GetString(object, tag, &dcminfo->AcquisitionTime);
  if (cond != DCM_NORMAL || strcmp(dcminfo->AcquisitionTime, "00000000")==0) {
    dcminfo->AcquisitionTime=(char *)calloc(256, sizeof(char));
    strcpy(dcminfo->AcquisitionTime, "UNKNOWN");
    cond2=cond;
#ifdef _DEBUG
    printf("WARNING: tag AcquisitionTime not found\n"); 
#endif  
  }
  else
    IsTagPresent[DCM_AcquisitionTime]=true;

  // slice thickness
  tag=DCM_MAKETAG(0x18, 0x50);
  cond=GetDoubleFromString(object, tag, &dcminfo->SliceThickness);
  if (cond != DCM_NORMAL || dcminfo->SliceThickness==0.0) {
    dcminfo->SliceThickness=0.0;
    cond2=cond;
#ifdef _DEBUG
    printf("WARNING: tag SliceThickness not found\n"); 
#endif  
  }
  else
    IsTagPresent[DCM_SliceThickness]=true;

  // image number
  tag=DCM_MAKETAG(0x20, 0x13);
  cond=GetUSFromString(object, tag, &dcminfo->ImageNumber);
  if ((cond != DCM_NORMAL || dcminfo->ImageNumber==0) && ImageNumber!=-1)
    {
      dcminfo->ImageNumber=ImageNumber;
      cond2=cond;
#ifdef _DEBUG
    printf("WARNING: tag ImageNumber not found\n"); 
#endif  
   }
  else
    IsTagPresent[DCM_ImageNumber]=true;

  // rows
  tag=DCM_MAKETAG(0x28, 0x10);
  cond=GetUSFromUS(object, tag, &dcminfo->Rows);
  if (cond != DCM_NORMAL) {
    dcminfo->Rows=0;
    cond2=cond;
#ifdef _DEBUG
    printf("WARNING: tag Rows not found\n"); 
#endif  
  }
  else
    IsTagPresent[DCM_Rows]=true;

  // columns
  tag=DCM_MAKETAG(0x28, 0x11);
  cond=GetUSFromUS(object, tag, &dcminfo->Columns);
  if (cond != DCM_NORMAL) {
    dcminfo->Columns=0;
    cond2=cond;
#ifdef _DEBUG
    printf("WARNING: tag Columns not found\n"); 
#endif  
  }
  else
    IsTagPresent[DCM_Columns]=true;

  // pixel spacing
  tag=DCM_MAKETAG(0x28, 0x30);
  cond=GetMultiDoubleFromString(object, tag, &tmp, 2);
  if (cond != DCM_NORMAL || tmp[0]==0.0 || tmp[1]==0.0) {
    dcminfo->xsize=0.0;
    dcminfo->ysize=0.0;
    cond2=cond;
#ifdef _DEBUG
    printf("WARNING: tag Pixel spacing not found\n"); 
#endif  
  }
  else {
    dcminfo->xsize=tmp[0];
    dcminfo->ysize=tmp[1];
    IsTagPresent[DCM_xsize]=true;
    IsTagPresent[DCM_ysize]=true;
  }

  // bits allocated
  tag=DCM_MAKETAG(0x28, 0x100);
  cond=GetUSFromUS(object, tag, &dcminfo->BitsAllocated);
  if (cond != DCM_NORMAL) {
    dcminfo->BitsAllocated=0;
    cond2=cond;
#ifdef _DEBUG
    printf("WARNING: tag BitsAllocated not found\n"); 
#endif  
  }
  else
    IsTagPresent[DCM_BitsAllocated]=true;

  // repetition time
  tag=DCM_MAKETAG(0x18, 0x80);
  cond=GetDoubleFromString(object, tag, &dcminfo->RepetitionTime);
  if (cond != DCM_NORMAL) 
    {
      dcminfo->RepetitionTime=0;
      cond2=cond;
#ifdef _DEBUG
    printf("WARNING: tag RepetitionTime not found\n"); 
#endif  
    }
  else
    IsTagPresent[DCM_RepetitionTime]=true;
  
  // echo time
  tag=DCM_MAKETAG(0x18, 0x81);
  cond=GetDoubleFromString(object, tag, &dcminfo->EchoTime);
  if (cond != DCM_NORMAL)
    {
      dcminfo->EchoTime=0;
      cond2=cond;
#ifdef _DEBUG
    printf("WARNING: tag EchoTime not found\n"); 
#endif  
    }
  else
    IsTagPresent[DCM_EchoTime]=true;

  // flip angle
  tag=DCM_MAKETAG(0x18, 0x1314);
  cond=GetDoubleFromString(object, tag, &dcminfo->FlipAngle);
  if (cond != DCM_NORMAL)
    {
      dcminfo->FlipAngle=0;
      cond2=cond;
#ifdef _DEBUG
    printf("WARNING: tag FlipAngle not found\n"); 
#endif  
    }
  else
    IsTagPresent[DCM_FlipAngle]=true;

  // inversion time
  tag=DCM_MAKETAG(0x18, 0x82);
  cond=GetDoubleFromString(object, tag, &dcminfo->InversionTime);
  if (cond != DCM_NORMAL)
    {
      dcminfo->InversionTime=0;
      cond2=cond;
#ifdef _DEBUG
    printf("WARNING: tag InversionTime not found\n"); 
#endif  
    }
  else
    IsTagPresent[DCM_InversionTime]=true;

  // echo number
  tag=DCM_MAKETAG(0x18, 0x86);
  cond=GetShortFromString(object, tag, &dcminfo->EchoNumber);
  if (cond != DCM_NORMAL)
    {
      dcminfo->EchoNumber=0;
      cond2=cond;
#ifdef _DEBUG
    printf("WARNING: tag EchoNumber not found\n"); 
#endif  
    }
  else
    IsTagPresent[DCM_EchoNumber]=true;

 // image position
  tag=DCM_MAKETAG(0x20, 0x32);
  cond=GetMultiDoubleFromString(object, tag, &tmp, 3);
  if (cond != DCM_NORMAL) 
    {
      for (i=0; i<3; i++)
  dcminfo->ImagePosition[i]=0;
      cond2=cond;
#ifdef _DEBUG
      printf("WARNING: tag image position not found\n"); 
#endif  
    }
  else 
    {
      IsTagPresent[DCM_ImagePosition]=true;
      for (i=0; i<3; i++)
  dcminfo->ImagePosition[i]=tmp[i];
    }
  
  // image orientation
  tag=DCM_MAKETAG(0x20, 0x37);
  cond=GetMultiDoubleFromString(object, tag, &tmp, 6);
  if (cond != DCM_NORMAL) {
    for (i=0; i<6; i++)
      dcminfo->ImageOrientation[i]=0;
    cond2=cond;
#ifdef _DEBUG
    printf("WARNING: tag image orientation not found\n"); 
#endif  
  }
  else {
    IsTagPresent[DCM_ImageOrientation]=true;
    for (i=0; i<6; i++)
      dcminfo->ImageOrientation[i]=tmp[i];
  }



  // pixel data
  if (ReadImage) {
    tag=DCM_MAKETAG(0x7FE0, 0x10);
    cond=GetPixelData(object, tag, &dcminfo->PixelData);
    if (cond != DCM_NORMAL)
      {
  dcminfo->PixelData=NULL;
  cond2=cond;
#ifdef _DEBUG
    printf("WARNING: tag pixel data not found\n"); 
#endif  
      }
  }

  cond=DCM_CloseObject(object);
  if (cond != DCM_NORMAL)
    cond2=cond;
  
  dcminfo->FieldOfView=dcminfo->xsize*dcminfo->Rows;

  free(tmp);
  free(itmp);
  return cond2;
}


/*******************************************************
   ReadDICOMImage
   Author: Sebastien Gicquel
   Date: 06/04/2001
   input: array of file names, number of elements
   output: pixels image (8 or 16 bits), and DICOM informations stored in first array element
*******************************************************/

void *ReadDICOMImage(int nfiles, DICOMInfo **aDicomInfo)
{
  int n, i, j, bitsAllocated, numberOfFrames, offset, nvox;
  DICOMInfo dcminfo;
  unsigned char *PixelData8, *v8=NULL;
  unsigned short *PixelData16, *v16=NULL;
  int npix;
  double firstPosition[3], lastPosition[3];
  unsigned short int min16=65535, max16=0;
  unsigned short int min8=255, max8=0;

  npix=aDicomInfo[0]->Columns*aDicomInfo[0]->Rows;
  numberOfFrames=nfiles;
  nvox=npix*numberOfFrames;
  aDicomInfo[0]->NumberOfFrames=numberOfFrames;
  bitsAllocated=aDicomInfo[0]->BitsAllocated;

  printf("reading DICOM image...\n");

  switch (bitsAllocated)
    {
    case 8:
      v8=(unsigned char *)calloc(nvox, sizeof(unsigned char));
      if (v8==NULL) {
  printf("ReadDICOMImage: allocation problem\n");
  exit(1);
      }
      break;
    case 16:
      v16=(unsigned short *)calloc(nvox, sizeof(unsigned short));
      if (v16==NULL) {
  printf("ReadDICOMImage: allocation problem\n");
  exit(1);
      }
      break;
    } // switch

  for (n=0; n<nfiles; n++) 
    {
      GetDICOMInfo(aDicomInfo[n]->FileName, &dcminfo, TRUE, n);
      
      if (n==0)
  for (i=0; i<3; i++)
    aDicomInfo[0]->FirstImagePosition[i]=dcminfo.ImagePosition[i];
      if (n==nfiles-1)
  for (i=0; i<3; i++)
    aDicomInfo[0]->LastImagePosition[i]=dcminfo.ImagePosition[i];
      
      offset=npix*n;
      
      switch (dcminfo.BitsAllocated) 
  {
  case 8:  
    PixelData8=(unsigned char*)dcminfo.PixelData;
    for (j=0; j<npix; j++)
      {
        v8[offset+j]=PixelData8[j];
        if (PixelData8[j]>max8)
    max8=PixelData8[j];
        else if (PixelData8[j]<min8)
    min8=PixelData8[j];
      }
    aDicomInfo[0]->min8=min8;
    aDicomInfo[0]->max8=max8;    
    free(PixelData8);
    break;
    
  case 16:  
    PixelData16=(unsigned short*)dcminfo.PixelData;
    for (j=0; j<npix; j++)
      {
        v16[offset+j]=PixelData16[j];
        if (PixelData16[j]>max16)
    max16=PixelData16[j];
        else if (PixelData16[j]<min16)
    min16=PixelData16[j];
      }
    aDicomInfo[0]->min16=min16;
    aDicomInfo[0]->max16=max16;
    free(PixelData16);
    break;
  }
    }

  for (i=0; i<3; i++)
    aDicomInfo[0]->ImagePosition[i]=(firstPosition[i]+lastPosition[i])/2.;

  switch (bitsAllocated)
    {
    case 8:
      return((void *)v8);
      break;
    case 16:
      return((void *)v16);
      break;
    default:
      return NULL;
    }
}


/*******************************************************
   SortFiles
   Author: Sebastien Gicquel
   Date: 06/04/2001
   input: array of file names, number of elements
   output: array of structures DICOMInfo sorted by study date and image number, number of studies encountered in the list of files
*******************************************************/

void SortFiles(char *fNames[], int nFiles, DICOMInfo ***ptrDicomArray, int *nStudies)
{
  int n, npermut;
  bool done;
  
  DICOMInfo **dicomArray, *storage;
  dicomArray=(DICOMInfo **)calloc(nFiles, sizeof(DICOMInfo *));
  *ptrDicomArray=dicomArray;
  
  for (n=0; n<nFiles; n++)
    {
      if ((dicomArray[n]=(DICOMInfo *)calloc(1, sizeof(DICOMInfo)))==NULL)
  {
    printf("DICOM conversion (SortFiles): can not allocate %d bytes\n", sizeof(DICOMInfo));
    exit(1);
  }
      GetDICOMInfo(fNames[n], dicomArray[n], FALSE, 1);
    }

  // sort by acquisition time, then image number
  done=false;
  while (!done)
    {
      npermut=0;
      for (n=0; n<nFiles-1; n++)
  if (strcmp(dicomArray[n]->AcquisitionTime, dicomArray[n+1]->AcquisitionTime)>0)  // 2nd time inferior to first
  {
      storage=dicomArray[n];
      dicomArray[n]=dicomArray[n+1];
      dicomArray[n+1]=storage;
      npermut++;
    }
      done=(npermut==0);      
    }  

  // sort by image number
  done=false;
  while (!done)
    {
      npermut=0;
      for (n=0; n<nFiles-1; n++)
  if (strcmp(dicomArray[n]->AcquisitionTime, dicomArray[n+1]->AcquisitionTime)==0
      && dicomArray[n]->ImageNumber>dicomArray[n+1]->ImageNumber)
    {
      storage=dicomArray[n];
      dicomArray[n]=dicomArray[n+1];
      dicomArray[n+1]=storage;
      npermut++;
    }
      done=(npermut==0);      
    }  

  // check studies number
  *nStudies=1;
  for (n=0; n<nFiles-1; n++)
    if (strcmp(dicomArray[n]->AcquisitionTime, dicomArray[n+1]->AcquisitionTime)!=0)
      (*nStudies)++;
}

/*******************************************************
   IsDICOM
   Author: Sebastien Gicquel
   Date: 06/04/2001
   input: DICOM file name
   output: true if file is DICOM part 10 (version 3.0) compliant
*******************************************************/

int IsDICOM(char *fname)
{

  unsigned long options=DCM_PART10FILE /*| DCM_ORDERLITTLEENDIAN |  DCM_ORDERBIGENDIAN | DCM_FORMATCONVERSION*/;

  DCM_OBJECT** object=(DCM_OBJECT**)calloc(1, sizeof(DCM_OBJECT*));
  CONDITION cond;

  cond=DCM_OpenFile(fname, options, object);

#ifdef TOTO
  switch (cond)
    {
    case DCM_NORMAL: printf("DCM_NORMAL\n"); break;
    case DCM_ILLEGALOPTION: printf("DCM_ILLEGALOPTION\n"); break;
    case DCM_OBJECTCREATEFAILED: printf("DCM_OBJECTCREATEFAILED\n"); break;
    case DCM_FILEOPENFAILED: printf("DCM_FILEOPENFAILED\n"); break;
    case DCM_FILEACCESSERROR: printf("DCM_FILEACCESSERROR\n"); break;
    case DCM_ELEMENTOUTOFORDER: printf("DCM_ELEMENTOUTOFORDER\n"); break;
    default: printf("don't know what's wrong...\n"); break;
    }
#endif

  DCM_CloseObject(object);

  return(cond == DCM_NORMAL);
}

/*******************************************************
   ScanDir
   Author: Sebastien Gicquel
   Date: 06/04/2001
   input: directory name
   output: array of files listed in directory, and number of files
*******************************************************/

#ifdef Solaris
/* added by kteich for solaris, since it doesn't have them by default. */
/* these funcs Copyright (c) Joerg-R. Hill, December 2000 */

int scandir(const char *dir, struct dirent ***namelist,
            int (*select)(const struct dirent *),
            int (*compar)(const struct dirent **, const struct dirent **))
{
  DIR *d;
  struct dirent *entry;
  register int i=0;
  size_t entrysize;

  if ((d=opendir(dir)) == NULL)
    return(-1);

  *namelist=NULL;
  while ((entry=readdir(d)) != NULL)
  {
    if (select == NULL || (select != NULL && (*select)(entry)))
    {
      *namelist=(struct dirent **)realloc((void *)(*namelist),
                                          (size_t)((i+1)*sizeof(struct dirent *)));
      if (*namelist == NULL) return(-1);
      entrysize=sizeof(struct dirent)-sizeof(entry->d_name)+strlen(entry->d_name)+1;
      (*namelist)[i]=(struct dirent *)malloc(entrysize);
      if ((*namelist)[i] == NULL) return(-1);
      memcpy((*namelist)[i], entry, entrysize);
      i++;
    }
  }
  if (closedir(d)) return(-1);
  if (i == 0) return(-1);
  if (compar != NULL)
    qsort((void *)(*namelist), (size_t)i, sizeof(struct dirent *), compar);
    
  return(i);
}

int alphasort(const struct dirent **a, const struct dirent **b)
{
  return(strcmp((*a)->d_name, (*b)->d_name));
}

#endif

int ScanDir(char *PathName, char ***FileNames, int *NumberOfFiles)
{
  char **pfn;
  struct dirent **NameList;
  int i, length, pathlength;

  pathlength=strlen(PathName);
  /* select all directory entries, and sort them by name */
  *NumberOfFiles = scandir(PathName, &NameList, 0, alphasort);

  if (*NumberOfFiles < 0)
    return -1;

  pfn = (char **)calloc(*NumberOfFiles, sizeof(char *));
  for (i=0; i<*NumberOfFiles; i++)
    {
      length=pathlength+strlen(NameList[i]->d_name)+1;
      pfn[i]=(char *)calloc(length, sizeof(char));
      sprintf(pfn[i], "%s%s", PathName, NameList[i]->d_name);
    }

  free(NameList);
  *FileNames=pfn;
  return 0;
}


/*******************************************************
   CleanFileNames
   Author: Sebastien Gicquel
   Date: 06/04/2001
   input: array of file names, number of elements
   output: array of file names that are DICOM 3.0 compliant 
*******************************************************/

int CleanFileNames(char **FileNames, int NumberOfDICOMFiles, char ***CleanedFileNames)
{
  
  char **pfn;
  int i, j, length;

  pfn = (char **)calloc(NumberOfDICOMFiles, sizeof(char *));
  for (i=0, j=0; j<NumberOfDICOMFiles; i++)
    {
      if (IsDICOM(FileNames[i]))
  {
    length=strlen(FileNames[i])+1;
    pfn[j]=(char *)calloc(length, sizeof(char));
    strcpy(pfn[j], FileNames[i]);
    j++;
  }
    }

  *CleanedFileNames=pfn;
  return DCM_NOERROR;
}

int round(double d)
{
  double c, f, ddown, dup; 

  c=ceil(d);
  f=floor(d);
  ddown=d-f;
  dup=c-d;
  if (ddown<dup)
    return (int)f;
  else
    return (int)c;
}

/*******************************************************
   RASFromOrientation
   Author: Sebastien Gicquel
   Date: 06/05/2001
   input: structure MRI, structure DICOMInfo
   output: fill in MRI structure RAS-related fields
   LIMITATION: only if patient was lying supine head first
*******************************************************/

int RASFromOrientation(MRI *mri, DICOMInfo *dcm)
{

/*
RAS coordinates inside scanner, if and only if patient is lying supine head first

-1  0  0
 0 -1  0 = M1
 0  0  1

Transformation matrix from DICOM volume to scanner is given by tag (20, 37) Image Orientation
= first row direction and first column direction in scanner referential
= M2

Transformation matrix from DICOM volume to MGH volume is given by

0 1  0
1 0  1 = M3
0 0 -1

RAS coordinates inside MRI structure are given by

M4 = M3.inv(M2).M1
*/

  MATRIX *RasScanner,
    *Scanner2dicom,
    *InvScanner2Dicom,
    /**Dicom2Mgh,*/
    *RasMgh,
    *VolumeCenterXyz,
    *VolumeCenterRas;
  int zdir=0, zor;
  double centerX, centerY, centerZ;

  RasScanner=MatrixAlloc(3, 3, 1);
  RasScanner->rptr[1][1]=-1.;
  RasScanner->rptr[2][2]=-1.;
  RasScanner->rptr[3][3]=1.;
  
  /*
  Dicom2Mgh=MatrixAlloc(3, 3, 1); 
  Dicom2Mgh->rptr[1][1]=1.;
  Dicom2Mgh->rptr[2][2]=1.;
  Dicom2Mgh->rptr[3][3]=1.;
  */

  Scanner2dicom=MatrixAlloc(3, 3, 1);
  Scanner2dicom->rptr[1][1]=round(dcm->ImageOrientation[0]);
  Scanner2dicom->rptr[2][1]=round(dcm->ImageOrientation[1]);
  Scanner2dicom->rptr[3][1]=round(dcm->ImageOrientation[2]);
  Scanner2dicom->rptr[1][2]=round(dcm->ImageOrientation[3]);
  Scanner2dicom->rptr[2][2]=round(dcm->ImageOrientation[4]);
  Scanner2dicom->rptr[3][2]=round(dcm->ImageOrientation[5]);

  if (round(dcm->ImageOrientation[0])!=0)
    {
      if (round(dcm->ImageOrientation[4])!=0)
  zdir=2;
      else
  zdir=1;
    }
  else if (round(dcm->ImageOrientation[1])!=0)
    {
      if (round(dcm->ImageOrientation[3])!=0)
  zdir=2;
      else
  zdir=0;
    }
  else if (round(dcm->ImageOrientation[2])!=0)
    {
      if (round(dcm->ImageOrientation[3])!=0)
  zdir=1;
      else
  zdir=0;
    }
  else
    printf("WARNING: DICOM orientation problems\n");

#ifdef _DEBUG
  printf("image orientation:\n");
  for (i=0; i<3; i++)
    printf("%g\t", dcm->ImageOrientation[i]);
  printf("\n");
  for (i=3; i<6; i++)
    printf("%g\t", dcm->ImageOrientation[i]);
  printf("\n");
  printf("first image position:\n");
  for (i=0; i<3; i++)
    printf("%g\t", dcm->FirstImagePosition[i]);
  printf("\n");
  printf("last image position:\n");
  for (i=0; i<3; i++)
    printf("%g\t", dcm->LastImagePosition[i]);
  printf("\n");
#endif


  if (dcm->LastImagePosition[zdir]>dcm->FirstImagePosition[zdir])
    zor=+1;
  else 
    zor=-1;
    
  Scanner2dicom->rptr[zdir+1][3]=zor;
#ifdef _DEBUG
  printf("matrix %s\n", "Scanner2dicom");
  MatrixPrint(stdout, Scanner2dicom);
#endif
  InvScanner2Dicom=MatrixInverse(Scanner2dicom, NULL);
#ifdef _DEBUG
  printf("matrix %s\n", "InvScanner2Dicom");
  MatrixPrint(stdout, InvScanner2Dicom);
#endif
  RasMgh=MatrixMultiply(InvScanner2Dicom, RasScanner, NULL);
#ifdef _DEBUG
  printf("matrix %s\n", "RasMg");
  MatrixPrint(stdout, RasMgh);
#endif
  
  mri->x_r=RasMgh->rptr[1][1];
  mri->y_r=RasMgh->rptr[2][1];
  mri->z_r=RasMgh->rptr[3][1];
  mri->x_a=RasMgh->rptr[1][2];
  mri->y_a=RasMgh->rptr[2][2];
  mri->z_a=RasMgh->rptr[3][2];
  mri->x_s=RasMgh->rptr[1][3];
  mri->y_s=RasMgh->rptr[2][3];
  mri->z_s=RasMgh->rptr[3][3];

  // DICOM tag (20, 32) (Image Position) gives the coordinates of the center of the first pixel transmitted, relatively to the scanner
 
  VolumeCenterXyz=MatrixAlloc(3, 1, 1);

  //  VolumeCenterXyz->rptr[1][1]=(dcm->FirstImagePosition[0]+dcm->LastImagePosition[0])/2.;
  centerX=(dcm->FirstImagePosition[0]+dcm->LastImagePosition[0])/2.;
  centerY=(dcm->FirstImagePosition[1]+dcm->LastImagePosition[1])/2.;
  centerZ=(dcm->FirstImagePosition[2]+dcm->LastImagePosition[2])/2.;

  switch (zdir)
    {
    case 0:  
      // slices along the x direction, sagittal
      // (x, y, z) scanner = (z, y, x) image
      centerY+=dcm->ysize*((double)dcm->Rows-0.5);
      centerZ+=dcm->xsize*((double)dcm->Columns-0.5);
      mri->slice_direction=MRI_SAGITTAL;
      break;
    case 1:  
      // slices along the y direction, coronal
      // (x, y, z) scanner = (x, z, y) image
      centerX+=dcm->xsize*((double)dcm->Rows-0.5);
      centerZ+=dcm->ysize*((double)dcm->Columns-0.5);
      mri->slice_direction=MRI_CORONAL;
      break;
    case 2:  
      // slices along the z direction, transaxial
      // (x, y, z) scanner = (x, y, z) image
      centerX+=dcm->xsize*((double)dcm->Rows-0.5);
      centerY+=dcm->ysize*((double)dcm->Columns-0.5);
      mri->slice_direction=MRI_HORIZONTAL;
      break;
    }

  VolumeCenterXyz->rptr[1][1]=centerX;
  VolumeCenterXyz->rptr[2][1]=centerY;
  VolumeCenterXyz->rptr[3][1]=centerZ;

  VolumeCenterRas=MatrixMultiply(RasScanner, VolumeCenterXyz, NULL);

  mri->c_r=VolumeCenterRas->rptr[1][1];
  mri->c_a=VolumeCenterRas->rptr[2][1];
  mri->c_s=VolumeCenterRas->rptr[3][1];
  
  mri->ras_good_flag=1;

  MatrixFree(&RasScanner);
  MatrixFree(&Scanner2dicom);
  MatrixFree(&InvScanner2Dicom);
  //  MatrixFree(&Dicom2Mgh);
  MatrixFree(&RasMgh);
  MatrixFree(&VolumeCenterXyz);
  MatrixFree(&VolumeCenterRas);
 
  return 0;
}

/*******************************************************
   DICOM16To8
   Author: Sebastien Gicquel
   Date: 06/06/2001
   input: array of 16-bit pixels, number of elements
   output: array of 8-bit pixels
*******************************************************/
unsigned char *DICOM16To8(unsigned short *v16, int nvox)
{
  unsigned char *v8;
  int i;
  double min16, max16, min8, max8, ratio;

  min8=0; 
  max8=255;
  
  v8=(unsigned char *)calloc(nvox, sizeof(unsigned char));
  if (v8==NULL)
    {
      printf("DICOMInfo2MRI: can't allocate %d bytes\n", nvox);
      exit(1);
    }
  
  for (i=0, min16=65535, max16=0; i<nvox; i++)  
    {
      if (v16[i]>max16)
  max16=(double)v16[i];
      if (v16[i]<min16)
  min16=(double)v16[i];
    }

  ratio = (max8-min8)/(max16-min16); 
  for (i=0; i<nvox; i++)  
      v8[i]=(unsigned short)((double)(v16[i])*ratio);
  
  return v8;
}
/*******************************************************
   DICOMInfo2MRI
   Author: Sebastien Gicquel
   Date: 06/05/2001
   input: structure DICOMInfo, pixel data 
   output: fill in MRI structure, including the image
*******************************************************/

int DICOMInfo2MRI(DICOMInfo *dcm, unsigned char *data, MRI *mri)
{
  int i, j, k, nvox;

  // fill in the fields
  strcpy(mri->fname, dcm->FileName);

  mri->width=dcm->Columns;
  mri->height=dcm->Rows;
  mri->depth=dcm->NumberOfFrames;
  nvox=dcm->Rows*dcm->Columns*dcm->NumberOfFrames;

  mri->fov=dcm->FieldOfView;
  mri->thick=dcm->SliceThickness;
  mri->xsize=dcm->xsize;
  mri->ysize=dcm->ysize;
  mri->zsize=dcm->SliceThickness;
  mri->nframes=1;

  mri->tr=dcm->RepetitionTime;
  mri->te=dcm->EchoTime;
  mri->ti=dcm->InversionTime;
  mri->flip_angle=dcm->FlipAngle;

  mri->imnr0=1;
  mri->imnr1=dcm->NumberOfFrames;

  RASFromOrientation(mri, dcm);
  
  for (k=0; k<dcm->NumberOfFrames; k++)
    for (j=0; j<dcm->Columns; j++)
      for (i=0; i<dcm->Rows; i++)
  mri->slices[k][j][i]=data[i+j*dcm->Rows+k*dcm->Rows*dcm->Columns];

  return 0;
}


/*******************************************************
   DICOMRead
   Author: Sebastien Gicquel
   Date: 06/04/2001
   input: directory name, boolean 
   output: MRI structure, including the image if input boolean is true
   FileName points to a DICOM filename. All DICOM files in same directory will be read
*******************************************************/

int DICOMRead(char *FileName, MRI **mri, int ReadImage)
{
  MRI *pmri;
  char **CleanedFileNames, **FileNames, *c, PathName[256];
  int i, NumberOfFiles, NumberOfDICOMFiles, nStudies, error;
  int length;
  DICOMInfo **aDicomInfo;
  unsigned char *v8=NULL;
  unsigned short *v16=NULL;

  for (i=0; i< NUMBEROFTAGS; i++)
    IsTagPresent[i]=false;

  // find pathname from filename
  c = strrchr(FileName, '/');
  if(c == NULL)
    PathName[0] = '\0';
  else
  {
    length = (int)(c - FileName);
    strncpy(PathName, FileName, length);
    PathName[length] = '/';
    PathName[length+1] = '\0';
  }

  // scan directory
  error=ScanDir(PathName, &FileNames, &NumberOfFiles);

  for (i=0, NumberOfDICOMFiles=0; i<NumberOfFiles; i++)
    {
      if (IsDICOM(FileNames[i]))
  {
    NumberOfDICOMFiles++;
  }
    }
  printf("%d DICOM 3.0 files in list\n", NumberOfDICOMFiles);

  // no DICOM files in directory
  if (NumberOfDICOMFiles==0)
    {
      printf("Sorry, no DICOM files found\n");
      exit(1);
    }

  // remove non DICOM file names
  CleanFileNames(FileNames, NumberOfDICOMFiles, &CleanedFileNames);
  //for (i=0; i<NumberOfFiles; i++)
  //free(FileNames[i]);
  //free(FileNames);
  
  // sort DICOM files by study date, then image number
  SortFiles(CleanedFileNames, NumberOfDICOMFiles, &aDicomInfo, &nStudies);
  if (nStudies>1)
    {
      printf("WARNING: %d different studies identified\n", nStudies);
      printf("Please, clen up directory\n");
      exit(1);
    }

  switch (aDicomInfo[0]->BitsAllocated)
    {
    case 8:
      v8=(unsigned char *)ReadDICOMImage(NumberOfDICOMFiles, aDicomInfo);
      break;
    case 16:
      v16=(unsigned short *)ReadDICOMImage(NumberOfDICOMFiles, aDicomInfo);
      v8=DICOM16To8(v16, aDicomInfo[0]->Columns*aDicomInfo[0]->Rows*aDicomInfo[0]->NumberOfFrames);
      free(v16);
    }
  
  // display only first DICOM header
  PrintDICOMInfo(aDicomInfo[0]);

  pmri=MRIallocSequence(aDicomInfo[0]->Columns, aDicomInfo[0]->Rows, aDicomInfo[0]->NumberOfFrames, MRI_UCHAR, 1);
  DICOMInfo2MRI(aDicomInfo[0], v8, pmri);
  *mri=pmri;


  return 0;

}
