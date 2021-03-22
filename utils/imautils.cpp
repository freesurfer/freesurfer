/**
 * @brief Siemens IMA file format utilities
 *
 */
/*
 * Original Author: Doug Greve
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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "fio.h"
#include "machine.h"

#include "imautils.h"

IMA_DICTIONARY_ENTRY ImaDictionary[NMAX_IMA_DICTIONARY] = {};
int nImaDictionary = 0, ImaDictionaryGood = 0;
const char *imaTypeString[6] = {"short", "int", "long", "float", "double", "string"};
int imaTypeSize[6] = {sizeof(short), sizeof(int), sizeof(long), sizeof(float), sizeof(double), sizeof(char)};

static int imaSetDictEntry(int nthEntry, const char *key, int offset, const char *typestring, int nitems);
static int imaGetKeyEntryNo(const char *key);

/*--------------------------------------------------------------------
  imaLoadVal() - loads a value of length nbytes*nitems from the file
  stream offset bytes from the beginning. If pval is non-null, the the
  data are copied into pval. Otherwise, a buffer of length
  nbytes*nitems is allocated and its pointer returned; in that case,
  the caller is responsible for dealloc. If the architecture is
  Intel0-based, the byte order is reversed based on nybtes.
  --------------------------------------------------------------------*/
void *imaLoadVal(FILE *imafp, int offset, int nbytes, int nitems, void *pval)
{
  int r;

  /* Go to the designated offset relative to the start of the file*/
  r = fseek(imafp, offset, SEEK_SET);
  if (r != 0) {
    printf("ERROR: imaLoadVal: could not seek to %d\n", offset);
    return (NULL);
  }
  /* Allocate memory if need be */
  if (pval == NULL) {
    pval = (void *)calloc(nbytes, nitems);
    if (pval == NULL) {
      printf("ERROR: imaLoadVal: could not alloc %d bytes\n", nbytes);
      return (NULL);
    }
  }
  /* Read in nitems of size nbytes */
  r = fread(pval, nbytes, nitems, imafp);
  if (r != nitems) {
    printf(
        "ERROR: imaLoadVal: could not read %d items of "
        "%d bytes from offset %d\n",
        nitems,
        nbytes,
        offset);
    return (NULL);
  }
  /* If this is a 486 arch, swap bytes because ima's are always
     created on a Sun */
  if (Arch486()) {
    r = ByteSwapBuf(pval, nitems, nbytes);
    if (r != 0) return (NULL);
  }

  return (pval);
}
/*--------------------------------------------------------------------
  MkImaDictionary() - Creates the IMA dictionary. Each entry is a
  key string that the describes the component of the IMA header,
  the number of bytes from the beginning of the file, the data type
  string, and, for strings, the number of characters in the string.
  The order is arbitrary.
  --------------------------------------------------------------------*/
void MkImaDictionary(void)
{
  // extern int nImaDictionary, ImaDictionaryGood;
  int n = 0;

  imaSetDictEntry(n++, "G08_Ide_StudyDate_Year", 0, "long", 1);
  imaSetDictEntry(n++, "G08_Ide_StudyDate_Month", 4, "long", 1);
  imaSetDictEntry(n++, "G08_Ide_StudyDate_Day", 8, "long", 1);
  imaSetDictEntry(n++, "G08_Ide_StudyTime_Hour", 36, "long", 1);
  imaSetDictEntry(n++, "G08_Ide_StudyTime_Minute", 40, "long", 1);
  imaSetDictEntry(n++, "G08_Ide_StudyTime_Second", 44, "long", 1);
  imaSetDictEntry(n++, "G08_Ide_Manufacturer", 96, "string", 27);
  imaSetDictEntry(n++, "G08_Ide_InstitutionID", 105, "string", 27);
  imaSetDictEntry(n++, "G08_Ide_ManufacturerModel", 281, "string", 27);
  imaSetDictEntry(n++, "G08_Ide_OperatorIdentification", 358, "string", 25);

  imaSetDictEntry(n++, "G10_Pat_PatientName", 768, "string", 27);
  imaSetDictEntry(n++, "G10_Pat_PatientBirthdate_Year", 808, "long", 1);
  imaSetDictEntry(n++, "G10_Pat_PatientBirthdate_Month", 812, "long", 1);
  imaSetDictEntry(n++, "G10_Pat_PatientBirthdate_Day", 816, "long", 1);
  imaSetDictEntry(n++, "G10_PatMod_PatientSex", 820, "int", 1);

  imaSetDictEntry(n++, "G18_Acq_SliceThickness", 1544, "double", 1);
  imaSetDictEntry(n++, "G18_Acq_EchoTime", 1568, "double", 1);
  imaSetDictEntry(n++, "G18_Acq_InversionTime", 1576, "double", 1);
  imaSetDictEntry(n++, "G18_Acq_NumberOfAverages", 1584, "long", 1);
  imaSetDictEntry(n++, "G18_Acq_EchoNumber", 1604, "long", 1);
  imaSetDictEntry(n++, "G18_Acq_DeviceSerialNumber", 1612, "string", 27);

  imaSetDictEntry(n++, "G19_Acq2_Mr_FlipAngle", 2112, "double", 1);
  imaSetDictEntry(n++, "G19_Acq3_Mr_MagneticFieldStrength", 2560, "double", 1);
  imaSetDictEntry(n++, "G19_Acq3_Mr_BaseRawMatrixSize", 2864, "int", 1);
  imaSetDictEntry(n++, "G19_Acq4_CM_ParameterFileName", 2944, "string", 65);
  imaSetDictEntry(n++, "G19_Acq4_CM_SequenceFileName", 3009, "string", 65);

  imaSetDictEntry(n++, "G20_Rel_Study", 3200, "long", 1);
  // imaSetDictEntry(n++,"G20_Rel_Acquisition",            3208,"long",1);
  imaSetDictEntry(n++, "G20_Rel_Image", 3212, "long", 1);
  // imaSetDictEntry(n++,"G20_Rel_AcquisitionsInSeries",   3292,"long",1);

  imaSetDictEntry(n++, "G21_Rel1_CM_FoV_Height", 3744, "double", 1);
  imaSetDictEntry(n++, "G21_Rel1_CM_FoV_Width", 3752, "double", 1);
  imaSetDictEntry(n++, "G21_Rel1_ImagePosition_Sag", 3768, "double", 1);
  imaSetDictEntry(n++, "G21_Rel1_ImagePosition_Cor", 3776, "double", 1);
  imaSetDictEntry(n++, "G21_Rel1_ImagePosition_Tra", 3784, "double", 1);
  imaSetDictEntry(n++, "G21_Rel1_CM_ImageNormal_Sag", 3792, "double", 1);
  imaSetDictEntry(n++, "G21_Rel1_CM_ImageNormal_Cor", 3800, "double", 1);
  imaSetDictEntry(n++, "G21_Rel1_CM_ImageNormal_Tra", 3808, "double", 1);
  imaSetDictEntry(n++, "G21_Rel1_CM_ImageColumn_Sag", 3832, "double", 1);
  imaSetDictEntry(n++, "G21_Rel1_CM_ImageColumn_Cor", 3840, "double", 1);
  imaSetDictEntry(n++, "G21_Rel1_CM_ImageColumn_Tra", 3848, "double", 1);
  imaSetDictEntry(n++, "G21_Rel1_CM_ImageRow_Sag", 3856, "double", 1);
  imaSetDictEntry(n++, "G21_Rel1_CM_ImageRow_Cor", 3864, "double", 1);
  imaSetDictEntry(n++, "G21_Rel1_CM_ImageRow_Tra", 3872, "double", 1);
  imaSetDictEntry(n++, "G21_Rel2_Mr_NumberOfSlicesNom", 4004, "long", 1);
  imaSetDictEntry(n++, "G21_Rel2_Mr_RepetitionTime", 4072, "double", 1);
  imaSetDictEntry(n++, "G21_Rel2_Mr_NumberOfEchoes", 4108, "long", 1);
  imaSetDictEntry(n++, "G21_Rel2_Mr_CurrentSliceDistanceFactor", 4136, "double", 1);

  imaSetDictEntry(n++, "G28_Pre_Rows", 4994, "short", 1);
  imaSetDictEntry(n++, "G28_Pre_Columns", 4996, "short", 1);
  imaSetDictEntry(n++, "G28_Pre_PixelSize_Row", 5000, "double", 1);
  imaSetDictEntry(n++, "G28_Pre_PixelSize_Column", 5008, "double", 1);
  imaSetDictEntry(n++, "G28_Pre_BitsAllocated", 5024, "short", 1);
  imaSetDictEntry(n++, "G28_Pre_BitsStored", 5026, "short", 1);
  imaSetDictEntry(n++, "G28_Pre_HighBit", 5028, "short", 1);

  imaSetDictEntry(n++, "G51_Txt_PatientSexAndAge", 5517, "string", 12);
  imaSetDictEntry(n++, "G51_Txt_PatientPosition", 5529, "string", 12);
  imaSetDictEntry(n++, "G51_Txt_ImageNumber", 5541, "string", 12);

  nImaDictionary = n;
  ImaDictionaryGood = 1;

  return;
}
/*--------------------------------------------------------------------
  --------------------------------------------------------------------*/
static int imaSetDictEntry(int nthEntry, const char *key, int offset, const char *typestring, int nitems)
{
  int type;

  if (nthEntry >= NMAX_IMA_DICTIONARY) {
    printf("ERROR: maximum number of ima Dictionary Entries exceeded\n");
    return (1);
  }

  ImaDictionary[nthEntry].key = key;
  ImaDictionary[nthEntry].offset = offset;
  ImaDictionary[nthEntry].typestring = typestring;
  type = imaTypeFromString(typestring);
  ImaDictionary[nthEntry].type = type;
  if (type < 0) {
    printf("ERROR: ima type is < 0\n");
    return (1);
  }
  ImaDictionary[nthEntry].typesize = imaTypeSize[type];

  if (type != IMA_TYPE_STRING)
    ImaDictionary[nthEntry].nitems = 1;
  else
    ImaDictionary[nthEntry].nitems = nitems;

  return (0);
}
/*---------------------------------------------------------------
  DumpImaDictionary() - dumps the dictionary to a file stream.
  This is just a list of the items in the dictionary, not the
  values from a file.
  ---------------------------------------------------------------*/
void DumpImaDictionary(FILE *fp)
{
  // extern IMA_DICTIONARY_ENTRY ImaDictionary[NMAX_IMA_DICTIONARY];
  // extern int nImaDictionary, ImaDictionaryGood;
  int n;
  IMA_DICTIONARY_ENTRY *ide;

  if (!ImaDictionaryGood) MkImaDictionary();

  for (n = 0; n < nImaDictionary; n++) {
    ide = &ImaDictionary[n];
    fprintf(
        fp, "%3d %-40s  %5d  %-7s %d  %3d\n", n, ide->key, ide->offset, ide->typestring, ide->typesize, ide->nitems);
  }
  return;
}
/*---------------------------------------------------------------
  DumpImaDictionaryVal() - dumps the dictionary to a file stream,
  including the corresponding values from a given ima file.
  ---------------------------------------------------------------*/
int DumpImaDictionaryVal(FILE *fp, const char *imafile)
{
  // extern IMA_DICTIONARY_ENTRY ImaDictionary[NMAX_IMA_DICTIONARY];
  // extern int nImaDictionary, ImaDictionaryGood;
  const char *key, *typestring;
  int n, offset, type, typesize, nitems;
  FILE *imafp;
  short sval;
  int ival;
  long lval;
  float fval;
  double dval;
  char string[1000];

  if (!ImaDictionaryGood) MkImaDictionary();

  imafp = fopen(imafile, "r");
  if (imafp == NULL) {
    printf("ERROR: cannot open %s\n", imafile);
    return (1);
  }

  for (n = 0; n < nImaDictionary; n++) {
    key = ImaDictionary[n].key;
    offset = ImaDictionary[n].offset;
    type = ImaDictionary[n].type;
    typesize = ImaDictionary[n].typesize;
    typestring = ImaDictionary[n].typestring;
    nitems = ImaDictionary[n].nitems;
    fprintf(fp, "%3d %-40s  %5d  %-7s  %3d   ", n, key, offset, typestring, nitems);

    switch (type) {
      case IMA_TYPE_SHORT:
        imaLoadVal(imafp, offset, typesize, nitems, &sval);
        fprintf(fp, "%d\n", (int)sval);
        break;
      case IMA_TYPE_INT:
        imaLoadVal(imafp, offset, typesize, nitems, &ival);
        fprintf(fp, "%d\n", ival);
        break;
      case IMA_TYPE_LONG:
        imaLoadVal(imafp, offset, typesize, nitems, &lval);
        fprintf(fp, "%d\n", (int)lval);
        break;
      case IMA_TYPE_FLOAT:
        imaLoadVal(imafp, offset, typesize, nitems, &fval);
        fprintf(fp, "%f\n", fval);
        break;
      case IMA_TYPE_DOUBLE:
        imaLoadVal(imafp, offset, typesize, nitems, &dval);
        fprintf(fp, "%f\n", (float)dval);
        break;
      case IMA_TYPE_STRING:
        imaLoadVal(imafp, offset, typesize, nitems, string);
        fprintf(fp, "%s\n", string);
        break;
    }
  }
  fclose(imafp);

  return (0);
}
/*--------------------------------------------------------------------
  imaPrintVal() - prints a value to a stream given a pointer to the
  value and the data type. The value is derefferenced based on the
  data type and printed. No spaces or newlines are printed.
  --------------------------------------------------------------------*/
int imaPrintVal(FILE *fp, int type, void *pval)
{
  switch (type) {
    case IMA_TYPE_SHORT:
      fprintf(fp, "%d", *((short *)pval));
      break;
    case IMA_TYPE_INT:
      fprintf(fp, "%d", *((int *)pval));
      break;
    case IMA_TYPE_LONG:
      fprintf(fp, "%ld", *((long *)pval));
      break;
    case IMA_TYPE_FLOAT:
      fprintf(fp, "%f", *((float *)pval));
      break;
    case IMA_TYPE_DOUBLE:
      fprintf(fp, "%f", (float)(*((double *)pval)));
      break;
    case IMA_TYPE_STRING:
      fprintf(fp, "%s", (char *)pval);
      break;
  }

  return (0);
}

/*--------------------------------------------------------------------
  imaLoadValFromKey() - loads a value from a file stream given the key
  in the IMA dictionary.
  --------------------------------------------------------------------*/
void *imaLoadValFromKey(FILE *imafp, const char *key, void *pval)
{
  // extern IMA_DICTIONARY_ENTRY ImaDictionary[NMAX_IMA_DICTIONARY];
  // extern int ImaDictionaryGood;
  int n, offset, typesize, nitems, nbytes;
  void *r;

  if (!ImaDictionaryGood) MkImaDictionary();

  n = imaGetKeyEntryNo(key);
  if (n < 0) return (NULL);

  offset = ImaDictionary[n].offset;
  typesize = ImaDictionary[n].typesize;
  nitems = ImaDictionary[n].nitems;
  nbytes = typesize * nitems;

  if (pval == NULL) {
    pval = (void *)calloc(nbytes, 1);
    if (pval == NULL) {
      printf("ERROR: imaLoadValbyKey: could not alloc %d bytes\n", nbytes);
      return (NULL);
    }
  }

  r = imaLoadVal(imafp, offset, typesize, nitems, pval);
  if (r == NULL) {
    free(pval);
    return (NULL);
  }

  return (pval);
}
/*--------------------------------------------------------------------
  imaTypeFromKey() - returns the interger code of the data type
  (IMA_TYPE_XXX) of an entry in the IMA dictionary given the key.
  --------------------------------------------------------------------*/
int imaTypeFromKey(const char *key)
{
  // extern IMA_DICTIONARY_ENTRY ImaDictionary[NMAX_IMA_DICTIONARY];
  // extern int  ImaDictionaryGood;
  int n;

  if (!ImaDictionaryGood) MkImaDictionary();

  n = imaGetKeyEntryNo(key);
  if (n < 0) return (-1);

  return (ImaDictionary[n].type);
}
/*--------------------------------------------------------------------
  imaGetKeyEntryNo() - return the entry number of the key.
  --------------------------------------------------------------------*/
static int imaGetKeyEntryNo(const char *key)
{
  // extern IMA_DICTIONARY_ENTRY ImaDictionary[NMAX_IMA_DICTIONARY];
  // extern int nImaDictionary, ImaDictionaryGood;
  int n;

  if (!ImaDictionaryGood) MkImaDictionary();

  for (n = 0; n < nImaDictionary; n++)
    if (strcmp(ImaDictionary[n].key, key) == 0) return (n);

  printf("ERROR: key %s not found\n", key);
  return (-1);
}
/*--------------------------------------------------------------------
  imaTypeFromString() - returns the interger type code (IMA_TYPE_XXX)
  given the string name of the type.
  --------------------------------------------------------------------*/
int imaTypeFromString(const char *typestring)
{
  int i;
  for (i = 0; i < 6; i++) {
    // printf("%d %s\n",i,imaTypeString[i]);
    if (strcmp(typestring, imaTypeString[i]) == 0) return (i);
  }
  printf("ERROR: no support for type %s\n", typestring);
  return (-1);
}
/*--------------------------------------------------------------------
  imaLoadFileInfo() - fills the IMAFILEINFO structure from the header
  of the given file. Determines whether the image is a mosaic. All
  times are converted to seconds. All angles to radians. XYZ values
  are converted to RAS. VolDim[2] for non-mosaics is set to -1. The
  distance between slices is computed from the slice thickness and
  the distance factor.
  --------------------------------------------------------------------*/
IMAFILEINFO *imaLoadFileInfo(const char *imafile)
{
  IMAFILEINFO *ifi;
  FILE *fp;
  int itmp, len, err;
  double dtmp, FoVHeight, FoVWidth;
  long ltmp, Year, Month, Day, Hour, Min, Sec;
  short stmp;
  char tmpstr[1000];
  int FirstImageNo;
  int nVolVoxs, nMosVoxs;
  char Separator;

  fp = fopen(imafile, "r");
  if (fp == NULL) {
    printf("ERROR: could not open %s\n", imafile);
    return (NULL);
  }

  ifi = (IMAFILEINFO *)calloc(sizeof(IMAFILEINFO), 1);

  err = imaParseName(imafile, &ifi->StudyNo, &ifi->SeriesNo, &ifi->ImageNo, &Separator);
  if (err) {
    free(ifi);
    return (NULL);
  }

  len = strlen(imafile);
  ifi->FileName = (char *)calloc(sizeof(char), len + 1);
  memmove(ifi->FileName, imafile, len);

  ifi->NFilesInSeries = imaCountFilesInSeries(imafile, &FirstImageNo);

  ifi->PatientName = (char *)imaLoadValFromKey(fp, "G10_Pat_PatientName", NULL);
  ifi->PulseSequence = (char *)imaLoadValFromKey(fp, "G19_Acq4_CM_SequenceFileName", NULL);

  imaLoadValFromKey(fp, "G10_Pat_PatientBirthdate_Year", &Year);
  imaLoadValFromKey(fp, "G10_Pat_PatientBirthdate_Month", &Month);
  imaLoadValFromKey(fp, "G10_Pat_PatientBirthdate_Day", &Day);
  sprintf(tmpstr, "%04ld%02ld%02ld", Year, Month, Day);
  len = strlen(tmpstr);
  ifi->PatientDOB = (char *)calloc(sizeof(char), len + 1);
  memmove(ifi->PatientDOB, tmpstr, len);

  imaLoadValFromKey(fp, "G10_PatMod_PatientSex", &itmp);
  if (itmp == 1)
    memmove(ifi->PatientGender, "female", 6);
  else
    memmove(ifi->PatientGender, "male", 4);

  imaLoadValFromKey(fp, "G08_Ide_StudyDate_Year", &Year);
  imaLoadValFromKey(fp, "G08_Ide_StudyDate_Month", &Month);
  imaLoadValFromKey(fp, "G08_Ide_StudyDate_Day", &Day);
  sprintf(tmpstr, "%04ld%02ld%02ld", Year, Month, Day);
  len = strlen(tmpstr);
  ifi->StudyDate = (char *)calloc(sizeof(char), len + 1);
  memmove(ifi->StudyDate, tmpstr, len);

  imaLoadValFromKey(fp, "G08_Ide_StudyTime_Hour", &Hour);
  imaLoadValFromKey(fp, "G08_Ide_StudyTime_Minute", &Min);
  imaLoadValFromKey(fp, "G08_Ide_StudyTime_Second", &Sec);
  sprintf(tmpstr, "%02ld%02ld%02ld", Hour, Min, Sec);
  len = strlen(tmpstr);
  ifi->StudyTime = (char *)calloc(sizeof(char), len + 1);
  memmove(ifi->StudyTime, tmpstr, len);

  imaLoadValFromKey(fp, "G19_Acq2_Mr_FlipAngle", &dtmp);
  ifi->FlipAngle = (float)M_PI * dtmp / 180.0; /* Convert to radians */
  imaLoadValFromKey(fp, "G18_Acq_EchoTime", &dtmp);
  ifi->EchoTime = (float)dtmp / 1000.0; /* Seconds */
  imaLoadValFromKey(fp, "G21_Rel2_Mr_RepetitionTime", &dtmp);
  ifi->RepetitionTime = (float)dtmp / 1000.0; /* Seconds */
  imaLoadValFromKey(fp, "G18_Acq_InversionTime", &dtmp);
  ifi->InversionTime = (float)dtmp / 1000.0; /* Seconds */

  /* Image Size (Mosaic Size for Functionals) */
  imaLoadValFromKey(fp, "G28_Pre_Rows", &stmp);
  ifi->NImageRows = (int)stmp;
  imaLoadValFromKey(fp, "G28_Pre_Columns", &stmp);
  ifi->NImageCols = (int)stmp;

  /* Number of cols and rows in volume */
  imaLoadValFromKey(fp, "G19_Acq3_Mr_BaseRawMatrixSize", &ifi->VolDim[0]);
  ifi->VolDim[1] = ifi->VolDim[0];

  /* Its a mosaic if the number of image cols is different from the
     number of volume cols */
  if (ifi->VolDim[0] != ifi->NImageCols)
    ifi->IsMosaic = 1;
  else
    ifi->IsMosaic = 0;

  /* If its not a mosaic, there is no way to tell the number of slices
     from the info in the file. If it is a mosaic, use the nominal */
  if (!ifi->IsMosaic)
    ifi->VolDim[2] = ifi->NFilesInSeries;
  else {
    imaLoadValFromKey(fp, "G21_Rel2_Mr_NumberOfSlicesNom", &ltmp);
    ifi->VolDim[2] = (int)ltmp;
  }

  /* Col and Row Spacing, should work for mosaics and non */
  imaLoadValFromKey(fp, "G21_Rel1_CM_FoV_Height", &dtmp);
  FoVHeight = (float)dtmp;
  imaLoadValFromKey(fp, "G21_Rel1_CM_FoV_Width", &dtmp);
  FoVWidth = (float)dtmp;
  ifi->VolRes[0] = FoVWidth / ifi->VolDim[0];  /* col res */
  ifi->VolRes[1] = FoVHeight / ifi->VolDim[1]; /* row res */

  /* Slice Spacing */
  imaLoadValFromKey(fp, "G18_Acq_SliceThickness", &dtmp);
  ifi->SliceThickness = (float)dtmp;
  imaLoadValFromKey(fp, "G21_Rel2_Mr_CurrentSliceDistanceFactor", &dtmp);
  ifi->DistanceFactor = (float)dtmp;
  if (ifi->DistanceFactor > 0 && ifi->DistanceFactor <= 1)
    ifi->VolRes[2] = ifi->SliceThickness * (1.0 + ifi->DistanceFactor);
  else
    ifi->VolRes[2] = ifi->SliceThickness;

  /* Center of the Image; convert to RAS */
  imaLoadValFromKey(fp, "G21_Rel1_ImagePosition_Sag", &dtmp);
  ifi->ImgCenter[0] = -dtmp; /* X,R */
  imaLoadValFromKey(fp, "G21_Rel1_ImagePosition_Cor", &dtmp);
  ifi->ImgCenter[1] = dtmp; /* Y,A */
  imaLoadValFromKey(fp, "G21_Rel1_ImagePosition_Tra", &dtmp);
  ifi->ImgCenter[2] = -dtmp; /* Z,S */

  /* Column Direction Cosines; convert to RAS */
  imaLoadValFromKey(fp, "G21_Rel1_CM_ImageColumn_Sag", &dtmp);
  ifi->Vc[0] = -dtmp; /* X,R */
  imaLoadValFromKey(fp, "G21_Rel1_CM_ImageColumn_Cor", &dtmp);
  ifi->Vc[1] = dtmp; /* Y,A */
  imaLoadValFromKey(fp, "G21_Rel1_CM_ImageColumn_Tra", &dtmp);
  ifi->Vc[1] = -dtmp; /* Z,S */

  /* Row Direction Cosines; convert to RAS */
  imaLoadValFromKey(fp, "G21_Rel1_CM_ImageRow_Sag", &dtmp);
  ifi->Vr[0] = -dtmp; /* X,R */
  imaLoadValFromKey(fp, "G21_Rel1_CM_ImageRow_Cor", &dtmp);
  ifi->Vr[1] = dtmp; /* Y,A */
  imaLoadValFromKey(fp, "G21_Rel1_CM_ImageRow_Tra", &dtmp);
  ifi->Vr[1] = -dtmp; /* Z,S */

  /* Slice Direction Cosines; convert to RAS */
  imaLoadValFromKey(fp, "G21_Rel1_CM_ImageNormal_Sag", &dtmp);
  ifi->Vs[0] = -dtmp; /* X,R */
  imaLoadValFromKey(fp, "G21_Rel1_CM_ImageNormal_Cor", &dtmp);
  ifi->Vs[1] = dtmp; /* Y,A */
  imaLoadValFromKey(fp, "G21_Rel1_CM_ImageNormal_Tra", &dtmp);
  ifi->Vs[1] = -dtmp; /* Z,S */

  if (!ifi->IsMosaic) {
    ifi->NFilesPerFrame = ifi->VolDim[2];
    ifi->NFrames = 1;
  }
  else {
    nVolVoxs = ifi->VolDim[0] * ifi->VolDim[1] * ifi->VolDim[2];
    nMosVoxs = ifi->NImageRows * ifi->NImageCols;
    ifi->NFilesPerFrame = (int)(ceil((float)nVolVoxs / nMosVoxs));
    /* Number of frames */
    imaLoadValFromKey(fp, "G18_Acq_NumberOfAverages", &ltmp);
    ifi->NFrames = (int)ltmp;
  }

  ifi->NFilesInSeriesExp = ifi->NFilesPerFrame * ifi->NFrames;

  if (ifi->NFilesInSeriesExp != ifi->NFilesInSeries)
    ifi->ErrorFlag = 1;
  else
    ifi->ErrorFlag = 0;

  fclose(fp);

  return (ifi);
}
/*--------------------------------------------------------------------*/
short *imaReadPixelData(IMAFILEINFO *ifi, short *PixelData)
{
  FILE *fp;
  int npixels, alloced = 0;
  int nread, r;

  fp = fopen(ifi->FileName, "r");
  if (fp == NULL) {
    printf("ERROR: cannot open %s\n", ifi->FileName);
    return (NULL);
  }

  npixels = ifi->NImageRows * ifi->NImageCols;

  if (PixelData == NULL) {
    PixelData = (short *)calloc(npixels, sizeof(short));
    if (PixelData == NULL) {
      printf("ERROR: could not alloc %d pixels\n", npixels);
      fclose(fp);
      return (NULL);
    }
    alloced = 1;
  }

  if (fseek(fp, 6144, SEEK_SET) == -1) {
    printf("ERROR: could seek to 6144 in %s\n", ifi->FileName);
    if (alloced) free(PixelData);
    fclose(fp);
    return (NULL);
  }

  nread = fread(PixelData, sizeof(short), npixels, fp);
  if (nread != npixels) {
    printf("ERROR: only read %d of %d pixels in  %s\n", nread, npixels, ifi->FileName);
    if (alloced) free(PixelData);
    fclose(fp);
    return (NULL);
  }
  fclose(fp);

  if (Arch486()) {
    r = ByteSwapBuf(PixelData, npixels, sizeof(short));
    if (r != 0) {
      if (alloced) free(PixelData);
      return (NULL);
    }
  }

  return (PixelData);
}
/*--------------------------------------------------------------------*/
int imaDumpFileInfo(FILE *fp, IMAFILEINFO *ifi)
{
  fprintf(fp, "FileName          %s\n", ifi->FileName);
  fprintf(fp, "PatientName       %s\n", ifi->PatientName);
  fprintf(fp, "PatientDOB        %s\n", ifi->PatientDOB);
  fprintf(fp, "PatientGender     %s\n", ifi->PatientGender);
  fprintf(fp, "StudyDate         %s\n", ifi->StudyDate);
  fprintf(fp, "StudyTime         %s\n", ifi->StudyTime);
  // fprintf(fp,"SeriesTime        %s\n",ifi->SeriesTime);
  // fprintf(fp,"AcqTime           %s\n",ifi->AcquisitionTime);
  fprintf(fp, "PulseSeq          %s\n", ifi->PulseSequence);
  // fprintf(fp,"Protocol          %s\n",ifi->ProtocolName);
  // fprintf(fp,"PhEncDir          %s\n",ifi->PhEncDir);
  // fprintf(fp,"EchoNo            %d\n",ifi->EchoNo);
  fprintf(fp, "FlipAngle         %g\n", ifi->FlipAngle);
  fprintf(fp, "EchoTime          %g\n", ifi->EchoTime);
  fprintf(fp, "InversionTime     %g\n", ifi->InversionTime);
  fprintf(fp, "RepetitionTime    %g\n", ifi->RepetitionTime);
  // fprintf(fp,"PhEncFOV          %g\n",ifi->PhEncFOV);
  // fprintf(fp,"ReadoutFOV        %g\n",ifi->ReadoutFOV);

  fprintf(fp, "SeriesNo          %d\n", ifi->SeriesNo);
  fprintf(fp, "NFilesInSeries    %d\n", ifi->NFilesInSeries);
  fprintf(fp, "ImageNo           %d\n", ifi->ImageNo);
  fprintf(fp, "NImageRows        %d\n", ifi->NImageRows);
  fprintf(fp, "NImageCols        %d\n", ifi->NImageCols);
  fprintf(fp, "NFrames           %d\n", ifi->NFrames);
  fprintf(fp, "IsMosaic          %d\n", ifi->IsMosaic);

  fprintf(fp, "ImgCenter  %8.4f %8.4f %8.4f \n", ifi->ImgCenter[0], ifi->ImgCenter[1], ifi->ImgCenter[2]);

  fprintf(fp, "VolRes     %8.4f %8.4f %8.4f \n", ifi->VolRes[0], ifi->VolRes[1], ifi->VolRes[2]);
  fprintf(fp, "VolDim     %3d      %3d       %3d \n", ifi->VolDim[0], ifi->VolDim[1], ifi->VolDim[2]);
  fprintf(fp, "Vc         %8.4f %8.4f %8.4f \n", ifi->Vc[0], ifi->Vc[1], ifi->Vc[2]);
  fprintf(fp, "Vr         %8.4f %8.4f %8.4f \n", ifi->Vr[0], ifi->Vr[1], ifi->Vr[2]);
  fprintf(fp, "Vs         %8.4f %8.4f %8.4f \n", ifi->Vs[0], ifi->Vs[1], ifi->Vs[2]);

  // fprintf(fp,"VolCenter  %8.4f %8.4f %8.4f \n",
  // ifi->VolCenter[0],ifi->VolCenter[1],ifi->VolCenter[2]);

  return (0);
}

/*--------------------------------------------------------------------
  imaParseName() - studyno-seriesno-imageno.ima
  --------------------------------------------------------------------*/
int imaParseName(const char *imafile, int *StudyNo, int *SeriesNo, int *ImageNo, char *Separator)
{
  char *imabase;
  int baselen, n, m;
  char tmpstr[500];

  imabase = fio_basename(imafile, NULL);
  baselen = strlen(imabase);

  // printf("%s\n",imafile);
  // printf("%s, %d\n",imabase,baselen);

  /* First try with Separator = '-' */
  *Separator = '-';

  /* Read the Series Number */
  memset(tmpstr, 0, 500);
  n = 0;
  m = 0;
  while (n < baselen && imabase[n] != *Separator) {
    // printf("%2d %c  %2d\n",n,imabase[n],m);
    tmpstr[m] = imabase[n];
    n++;
    m++;
  }
  if (n == baselen) {
    /* '-' did not work, try with '_' */
    *Separator = '_';
    n = 0;
    m = 0;
    while (n < baselen && imabase[n] != *Separator) {
      // printf("%2d %c  %2d\n",n,imabase[n],m);
      tmpstr[m] = imabase[n];
      n++;
      m++;
    }
    if (n == baselen) {
      printf("ERROR: could not parse %s\n", imafile);
      return (1);
    }
  }
  sscanf(tmpstr, "%d", StudyNo);
  // printf("StudyNo %d\n",*StudyNo);

  /* Read the Series (Run) Number */
  memset(tmpstr, 0, 500);
  n++;
  m = 0;
  while (n < baselen && imabase[n] != *Separator) {
    tmpstr[m] = imabase[n];
    n++;
    m++;
  }
  if (n == baselen) {
    printf("ERROR: could not parse %s\n", imafile);
    return (1);
  }
  sscanf(tmpstr, "%d", SeriesNo);
  // printf("Series %d\n",*SeriesNo);

  /* Read the Image Number */
  memset(tmpstr, 0, 500);
  n++;
  m = 0;
  while (n < baselen && imabase[n] != '.') {
    tmpstr[m] = imabase[n];
    n++;
    m++;
  }
  if (n == baselen) {
    printf("ERROR: could not parse %s\n", imafile);
    return (1);
  }
  sscanf(tmpstr, "%d", ImageNo);
  // printf("Image %d\n",*ImageNo);

  free(imabase);

  return (0);
}
/*--------------------------------------------------------------------
  imaHasIMAExtension() - returns 1 if the given file name ends in a
  .ima (or .IMA) extension.
  --------------------------------------------------------------------*/
int imaHasIMAExtension(const char *filename)
{
  char *ext;

  ext = fio_extension(filename);
  if (ext == NULL) return (0);

  if (!strcasecmp(ext, "ima")) return (1);

  return (0);
}
/*--------------------------------------------------------------------
  imaIsSiemensIMA() - returns 1 if the given file is a Siemens IMA
  file. Returns 0 if it is not a Siemens IMA file or if the file does
  not exist.
  --------------------------------------------------------------------*/
int imaIsSiemensIMA(const char *imafile)
{
  FILE *fp;
  char *Manufacturer;

  fp = fopen(imafile, "r");
  if (fp == NULL) {
    // printf("ERROR: could not open %s\n",imafile);
    return (0);
  }

  Manufacturer = (char *)imaLoadValFromKey(fp, "G08_Ide_Manufacturer", NULL);

  if (Manufacturer == NULL) return (0);
  if (!strcasecmp(Manufacturer, "SIEMENS")) {
    free(Manufacturer);
    return (1);
  }

  free(Manufacturer);
  return (0);
}
/*--------------------------------------------------------------------
  imaCountFilesInSeries() - Assuming that the IMA file name is
  studyno-seriesno-imageno.ima, this will look for files
  down from imageno to determine the image number of the first
  file in the series, and up to determine the last. The number
  in the series is then last-first+1.
  --------------------------------------------------------------------*/
int imaCountFilesInSeries(const char *imafile, int *FirstImageNo)
{
  int StudyNo, SeriesNo, ImageNo, LastImageNo;
  char *imadir = NULL, *imaext;
  int err, n;
  FILE *fp;
  char filename[1000];
  char Separator;

  fp = fopen(imafile, "r");
  if (fp == NULL) {
    printf("ERROR: could not open %s\n", imafile);
    return (0);
  }

  err = imaParseName(imafile, &StudyNo, &SeriesNo, &ImageNo, &Separator);
  if (err) return (-1);

  imadir = fio_dirname(imafile);
  imaext = fio_extension(imafile);

  /* Count down */
  n = ImageNo;
  while (1) {
    n--;
    if (Separator == '-')
      sprintf(filename, "%s/%d-%d-%d.%s", imadir, StudyNo, SeriesNo, n, imaext);
    else
      sprintf(filename, "%s/%d_%d_%d.%s", imadir, StudyNo, SeriesNo, n, imaext);
    // printf("%s\n",filename);
    fp = fopen(filename, "r");
    if (fp == NULL) break;
    fclose(fp);
  }
  *FirstImageNo = n + 1;

  /* Count up */
  n = ImageNo;
  while (1) {
    n++;
    if (Separator == '-')
      sprintf(filename, "%s/%d-%d-%d.%s", imadir, StudyNo, SeriesNo, n, imaext);
    else
      sprintf(filename, "%s/%d_%d_%d.%s", imadir, StudyNo, SeriesNo, n, imaext);
    // printf("%s\n",filename);
    fp = fopen(filename, "r");
    if (fp == NULL) break;
    fclose(fp);
  }
  LastImageNo = n - 1;

  free(imadir);
  free(imaext);

  return (LastImageNo - *FirstImageNo + 1);
}

/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/
#if 0
int imaAddElement(IMAFILEINFO *ifi, IMAELEMENT *e)
{
  IMAELEMENT **tmp;

  ifi->nelements ++;

  tmp = (IMAELEMENT **)realloc(ifi->e, ifi->nelements * sizeof(IMAELEMENT*) );
  if (tmp == NULL)
  {
    printf("ERROR: imaAddElement: could not alloc %d\n",ifi->nelements);
    ifi->nelements --;
    return(0);
  }
  ifi->e = tmp;
  ifi->e[ifi->nelements-1] = imaCopyElement(e);

  //printf("nelements = %d\n",ifi->nelements);

  return(0);
}
/*--------------------------------------------------------------------*/
int imaDumpFileInfo(FILE *fp, IMAFILEINFO *ifi)
{
  int n;

  fprintf(fp,"filename %s\n",ifi->filename);
  fprintf(fp,"nelements %d\n",ifi->nelements);

  for (n=0; n < ifi->nelements; n++)
  {
    printf("%3d ",n);
    imaDumpElement(fp,ifi->e[n]);
  }

  return(0);
}
/*---------------------------------------------------------*/
IMAFILEINFO *imaLoadDefaultFileInfo(const char *imafile)
{
  FILE *fp;
  IMAFILEINFO *ifi;
  int err,n;

  fp = fopen(imafile,"r");
  if (fp == NULL)
  {
    printf("ERROR: could not open %s for reading\n",imafile);
    return(NULL);
  }

  ifi = imaDefaultFileInfo();

  n = strlen(imafile);
  ifi->filename = (char *) calloc(n+1,sizeof(char));
  memmove(ifi->filename,imafile,n);

  for (n=0; n < ifi->nelements; n++)
  {
    err = imaLoadElementVal(fp, ifi->e[n]);
    if (err) return(NULL);
  }

  return(ifi);
}
/*--------------------------------------------------------------------*/
IMAFILEINFO *imaDefaultFileInfo(void)
{
  IMAFILEINFO *ifi;
  IMAELEMENT *e;

  ifi = (IMAFILEINFO *) calloc(sizeof(IMAFILEINFO),1);
  ifi->filename = "nofile";

  e = imaMakeElement("nimgrows",4994,IMA_TYPE_SHORT,1);
  imaAddElement(ifi,e);
  imaFreeElement(&e);

  e = imaMakeElement("nimgcols",4996,IMA_TYPE_SHORT,1);
  imaAddElement(ifi,e);
  imaFreeElement(&e);

  e = imaMakeElement("rawsize",2864,IMA_TYPE_INT,1);
  imaAddElement(ifi,e);
  imaFreeElement(&e);

  e = imaMakeElement("G21_Rel2_Mr_NumberOfSlicesNom",4004,IMA_TYPE_LONG,1);
  imaAddElement(ifi,e);
  imaFreeElement(&e);
  e = imaMakeElement("G21_Rel2_Mr_NumberOfSlicesCur",4008,IMA_TYPE_LONG,1);
  imaAddElement(ifi,e);
  imaFreeElement(&e);
  e = imaMakeElement("G21_Rel2_Mr_CurrentSliceNumber",4012,IMA_TYPE_LONG,1);
  imaAddElement(ifi,e);
  imaFreeElement(&e);

  e = imaMakeElement("G18_Acq_InversionTime",1584,IMA_TYPE_INT,1);
  imaAddElement(ifi,e);
  imaFreeElement(&e);

  e = imaMakeElement("image_no",5541,IMA_TYPE_STRING,12);
  imaAddElement(ifi,e);
  imaFreeElement(&e);

  e = imaMakeElement("colres",5008,IMA_TYPE_DOUBLE,1);
  imaAddElement(ifi,e);
  imaFreeElement(&e);

  e = imaMakeElement("rowres",5000,IMA_TYPE_DOUBLE,1);
  imaAddElement(ifi,e);
  imaFreeElement(&e);

  e = imaMakeElement("G18_Acq_SliceThickness",1544,IMA_TYPE_DOUBLE,1);
  imaAddElement(ifi,e);
  imaFreeElement(&e);

  e = imaMakeElement("distance_factor",4136,IMA_TYPE_DOUBLE,1);
  imaAddElement(ifi,e);
  imaFreeElement(&e);

  e = imaMakeElement("repetition_time",1560,IMA_TYPE_DOUBLE,1);
  imaAddElement(ifi,e);
  imaFreeElement(&e);

  e = imaMakeElement("echo_time",1568,IMA_TYPE_DOUBLE,1);
  imaAddElement(ifi,e);
  imaFreeElement(&e);

  e = imaMakeElement("G21_Rel2_Mr_NumberOfEchoes",4108,IMA_TYPE_DOUBLE,1);
  imaAddElement(ifi,e);
  imaFreeElement(&e);

  e = imaMakeElement("inversion_time",1576,IMA_TYPE_DOUBLE,1);
  imaAddElement(ifi,e);
  imaFreeElement(&e);

  e = imaMakeElement("flip_angle",2112,IMA_TYPE_DOUBLE,1);
  imaAddElement(ifi,e);
  imaFreeElement(&e);

  e = imaMakeElement("G21_Rel1_ImagePosition_Sag",3768,IMA_TYPE_DOUBLE,1);
  imaAddElement(ifi,e);
  imaFreeElement(&e);
  e = imaMakeElement("G21_Rel1_ImagePosition_Cor",3776,IMA_TYPE_DOUBLE,1);
  imaAddElement(ifi,e);
  imaFreeElement(&e);
  e = imaMakeElement("G21_Rel1_ImagePosition_Tra",3784,IMA_TYPE_DOUBLE,1);
  imaAddElement(ifi,e);
  imaFreeElement(&e);

  e = imaMakeElement("col_dircos_x",3832,IMA_TYPE_DOUBLE,1);
  imaAddElement(ifi,e);
  imaFreeElement(&e);
  e = imaMakeElement("col_dircos_y",3840,IMA_TYPE_DOUBLE,1);
  imaAddElement(ifi,e);
  imaFreeElement(&e);
  e = imaMakeElement("col_dircos_z",3848,IMA_TYPE_DOUBLE,1);
  imaAddElement(ifi,e);
  imaFreeElement(&e);

  e = imaMakeElement("row_dircos_x",3856,IMA_TYPE_DOUBLE,1);
  imaAddElement(ifi,e);
  imaFreeElement(&e);
  e = imaMakeElement("row_dircos_y",3864,IMA_TYPE_DOUBLE,1);
  imaAddElement(ifi,e);
  imaFreeElement(&e);
  e = imaMakeElement("row_dircos_z",3872,IMA_TYPE_DOUBLE,1);
  imaAddElement(ifi,e);
  imaFreeElement(&e);

  e = imaMakeElement("slice_dircos_x",3792,IMA_TYPE_DOUBLE,1);
  imaAddElement(ifi,e);
  imaFreeElement(&e);
  e = imaMakeElement("slice_dircos_y",3800,IMA_TYPE_DOUBLE,1);
  imaAddElement(ifi,e);
  imaFreeElement(&e);
  e = imaMakeElement("slice_dircos_z",3808,IMA_TYPE_DOUBLE,1);
  imaAddElement(ifi,e);
  imaFreeElement(&e);

  //e = imaMakeElement("G21_Rel2_Mr_PhaseEncodingVectorTra",4160,IMA_TYPE_FLOAT,1);
  //imaAddElement(ifi,e); imaFreeElement(&e);
  //e = imaMakeElement("G21_Rel2_Mr_ReadoutVectorTra",4184,IMA_TYPE_FLOAT,1);
  //imaAddElement(ifi,e); imaFreeElement(&e);



  e = imaMakeElement("G21_Rel1_CM_FoV_Width",3752,IMA_TYPE_DOUBLE,1);
  imaAddElement(ifi,e);
  imaFreeElement(&e);

  e = imaMakeElement("G21_Rel1_CM_FoV_Height",3744,IMA_TYPE_DOUBLE,1);
  imaAddElement(ifi,e);
  imaFreeElement(&e);


  //e = imaMakeElement("G51_Txt_FieldOfView",5838,IMA_TYPE_STRING,1);
  //imaAddElement(ifi,e); imaFreeElement(&e);

  //e = imaMakeElement("Rel1_CM_ImagePosition_Sag",3768,IMA_TYPE_DOUBLE,1);
  //imaAddElement(ifi,e); imaFreeElement(&e);

  //e = imaMakeElement("Pre_ImageDimension",4992,IMA_TYPE_INT,1);
  //imaAddElement(ifi,e); imaFreeElement(&e);


  //imaDumpFileInfo(stdout,ifi);

  return(ifi);
}
/*--------------------------------------------------------------------*/
int imaDumpElement(FILE *fp, IMAELEMENT *e)
{
  fprintf(fp,"%-28s %5d %-6s %2d : ",
          e->descr, e->offset,
          imaTypeString[e->type], e->nbytes*e->nitems);

  if (e->valgood) imaPrintElementValue(fp, e);
  fprintf(fp,"\n");

  return(0);
}
/*--------------------------------------------------------------------*/
int imaPrintElementValue(FILE *fp, IMAELEMENT *e)
{
  switch (e->type)
  {
  case IMA_TYPE_SHORT:
    fprintf(fp,"%d",e->sval);
    break;
  case IMA_TYPE_INT:
    fprintf(fp,"%d",e->ival);
    break;
  case IMA_TYPE_LONG:
    fprintf(fp,"%ld",e->lval);
    break;
  case IMA_TYPE_FLOAT:
    fprintf(fp,"%f",e->fval);
    break;
  case IMA_TYPE_DOUBLE:
    fprintf(fp,"%f",(float)e->dval);
    break;
  case IMA_TYPE_STRING:
    fprintf(fp,"%s",e->string);
    break;
  }
  return(0);
}

/*--------------------------------------------------------------------*/
int imaFreeElement(IMAELEMENT **ppe)
{
  IMAELEMENT *pe;
  pe = *ppe;

  free(pe->descr);
  if (pe->type == IMA_TYPE_STRING && pe->string != NULL)
    free(pe->string);
  free(*ppe);

  return(0);
}
/*--------------------------------------------------------------------*/
IMAELEMENT *imaMakeElement(const char *descr, int offset, int type, int nitems)
{
  IMAELEMENT *e;
  int descrlen;

  if (type < 0 || type > 5)
  {
    printf("ERROR: type = %d, must be between 0 and 5, inclusive\n",type);
    return(NULL);
  }
  if (descr == NULL)
  {
    printf("ERROR: description is NULL\n");
    return(NULL);
  }

  e = (IMAELEMENT *) calloc(1, sizeof(IMAELEMENT));

  e->offset = offset;
  e->type   = type;
  e->nbytes = imaTypeSize[type];
  e->valgood = 0;

  if (type != IMA_TYPE_STRING) e->nitems = 1;
  else                        e->nitems = nitems;

#if 0
  if (val != NULL)
  {
    if (type != IMA_TYPE_STRING)
    {
      switch (type)
      {
      case IMA_TYPE_SHORT:
        e->sval = *( (short *)  val);
        break;
      case IMA_TYPE_INT:
        e->ival = *( (int *)    val);
        break;
      case IMA_TYPE_LONG:
        e->lval = *( (long *)   val);
        break;
      case IMA_TYPE_FLOAT:
        e->fval = *( (float *)  val);
        break;
      case IMA_TYPE_DOUBLE:
        e->dval = *( (double *) val);
        break;
      }
    }
    else
    {
      e->string = (char *) calloc(e->nbytes+1,sizeof(char));
      memmove(e->string,val,e->nbytes);
    }
  }
#endif

  descrlen = strlen(descr);
  e->descr = (char *) calloc(descrlen+1,sizeof(char));
  memmove(e->descr,descr,descrlen);

  //printf("element descr: %s\n",e->descr);

  return(e);
}
/*--------------------------------------------------------------------*/
IMAELEMENT *imaCopyElement(IMAELEMENT *esrc)
{
  IMAELEMENT *edest;
  int descrlen;

  edest = (IMAELEMENT *) calloc(1, sizeof(IMAELEMENT));

  edest->offset  = esrc->offset;
  edest->type    = esrc->type;
  edest->nbytes  = esrc->nbytes;
  edest->nitems  = esrc->nitems;
  edest->valgood = esrc->valgood;
  edest->sval = esrc->sval;
  edest->ival = esrc->ival;
  edest->lval = esrc->lval;
  edest->fval = esrc->fval;
  edest->dval = esrc->dval;

  if (esrc->type == IMA_TYPE_STRING && esrc->nitems > 0 &&
      esrc->string != NULL)
  {
    edest->string = (char *) calloc(edest->nitems+1,sizeof(char));
    memmove(edest->string,esrc->string,edest->nitems);
  }
  descrlen = strlen(esrc->descr);
  edest->descr = (char *) calloc(descrlen+1,sizeof(char));
  memmove(edest->descr,esrc->descr,descrlen);

  return(edest);
}
/*--------------------------------------------------------------------*/
IMAELEMENT *imaLoadElement(FILE *imafp, int offset, int type,const char *descr,
                           int nitems)
{
  IMAELEMENT *e;
  int err;

  e = imaMakeElement(descr,offset,type,nitems);
  if (e == NULL) return(NULL);

  err = imaLoadElementVal(imafp, e);
  if (err)
  {
    imaFreeElement(&e);
    return(NULL);
  }
  return(e);
}

/*--------------------------------------------------------------------*/
int imaLoadElementVal(FILE *imafp, IMAELEMENT *e)
{
  void *v;
  int of, nb, ni;
  of = e->offset;
  nb = e->nbytes;
  ni = e->nitems;

  v = NULL;

  switch (e->type)
  {
  case IMA_TYPE_SHORT:
    v = imaLoadVal(imafp,of,nb,ni,(void *)(&e->sval));
    break;
  case IMA_TYPE_INT:
    v = imaLoadVal(imafp,of,nb,ni,(void *)(&e->ival));
    break;
  case IMA_TYPE_LONG:
    v = imaLoadVal(imafp,of,nb,ni,(void *)(&e->lval));
    break;
  case IMA_TYPE_FLOAT:
    v = imaLoadVal(imafp,of,nb,ni,(void *)(&e->fval));
    break;
  case IMA_TYPE_DOUBLE:
    v = imaLoadVal(imafp,of,nb,ni,(void *)(&e->dval));
    break;
  case IMA_TYPE_STRING:
    v = imaLoadVal(imafp,of,nb,ni,NULL);
    e->string = (char *) v;
    break;
  }

  if (v==NULL) return(1);

  e->valgood = 1;
  return(0);
}
#endif

#if 0
/*From:  ~rhoge/sde/minc/mgh/conversion/siemens_transfer/siemens_header_table.h*/
Siemens_header_entry Siemens_header_table[] = {
  {0x0008, 0x0020, &Siemens_header.G08.Ide.StudyDate, create_ds_date_t_element, 1},
  {0x0008, 0x0022, &Siemens_header.G08.Ide.AcquisitionDate, create_ds_date_t_element, 1},
  {0x0008, 0x0023, &Siemens_header.G08.Ide.ImageDate, create_ds_date_t_element, 1},
  {0x0008, 0x0030, &Siemens_header.G08.Ide.StudyTime, create_ds_time_t_element, 1},
  {0x0008, 0x0032, &Siemens_header.G08.Ide.AcquisitionTime, create_ds_time_t_element, 1},
  {0x0008, 0x0033, &Siemens_header.G08.Ide.ImageTime, create_ds_time_t_element, 1},
  {0x0008, 0x0041, &Siemens_header.G08.Ide.DataSetSubtype, create_data_set_subtype_t_element, 1},
  {0x0008, 0x0060, &Siemens_header.G08.Ide.Modality, create_modality_t_element, 1},
  {0x0008, 0x0070, &Siemens_header.G08.Ide.Manufacturer, create_char_element, LENGTH_MANUFACTURER + 1},
  {0x0008, 0x0080, &Siemens_header.G08.Ide.InstitutionID, create_char_element, LENGTH_LABEL + 1},
  {0x0008, 0x0090, &Siemens_header.G08.Ide.ReferringPhysician, create_char_element, LENGTH_LABEL + 1},
  {0x0008, 0x1010, &Siemens_header.G08.Ide.StationID, create_char_element, LENGTH_LABEL + 1},
  {0x0008, 0x1080, &Siemens_header.G08.Ide.AdmittingDiagnosis, create_char_element, LENGTH_DIAGNOSIS + 1},
  {0x0008, 0x1090, &Siemens_header.G08.Ide.ManufacturerModel, create_char_element, LENGTH_LABEL + 1},
  {0x0009, 0x1041, &Siemens_header.G09.Ide.DataObjectSubtype, create_data_object_subtype_t_element, 1},
  {0x0009, 0x1210, &Siemens_header.G09.Ide.StorageMode, create_storage_mode_t_element, 1},
  {0x0009, 0x1226, &Siemens_header.G09.Ide.LastMoveDate, create_ds_date_t_element, 1},
  {0x0009, 0x1227, &Siemens_header.G09.Ide.LastMoveTime, create_ds_time_t_element, 1},
  {0x0009, 0x1316, &Siemens_header.G09.Ide.CPUIdentificationLabel, create_char_element, LENGTH_LABEL + 1},
  {0x0009, 0x1320, &Siemens_header.G09.Ide.HeaderVersion, create_char_element, LENGTH_HEADER_VERSION + 1},
  {0x0010, 0x0010, &Siemens_header.G10.Pat.PatientName, create_char_element, LENGTH_LABEL + 1},
  {0x0010, 0x0020, &Siemens_header.G10.Pat.PatientId, create_char_element, LENGTH_PATIENT_ID + 1},
  {0x0010, 0x0030, &Siemens_header.G10.Pat.PatientBirthdate, create_ds_date_t_element, 1},
  {0x0010, 0x0040, &Siemens_header.G10.Pat.PatientSex, create_sex_t_element, 1},
  {0x0010, 0x1010, &Siemens_header.G10.Pat.PatientAge, create_char_element, LENGTH_AGE + 1},
  {0x0010, 0x1030, &Siemens_header.G10.Pat.PatientWeight, create_long_element, 1},
  {0x0011, 0x1110, &Siemens_header.G11.Pat.RegistrationDate, create_ds_date_t_element, 1},
  {0x0011, 0x1111, &Siemens_header.G11.Pat.RegistrationTime, create_ds_time_t_element, 1},
  {0x0011, 0x1123, &Siemens_header.G11.Pat.UsedPatientWeight, create_long_element, 1},
  {0x0018, 0x0010, &Siemens_header.G18.Acq.Contrast, create_contrast_t_element, 1},
  {0x0018, 0x0050, &Siemens_header.G18.Acq.SliceThickness, create_double_element, 1},
  {0x0018, 0x0080, &Siemens_header.G18.Acq.RepetitionTime, create_double_element, 1},
  {0x0018, 0x0081, &Siemens_header.G18.Acq.EchoTime, create_double_element, 1},
  {0x0018, 0x0083, &Siemens_header.G18.Acq.NumberOfAverages, create_long_element, 1},
  {0x0018, 0x0084, &Siemens_header.G18.Acq.ImagingFrequency, create_double_element, 1},
  {0x0018, 0x0085, &Siemens_header.G18.Acq.ImagedNucleus, create_char_element, LENGTH_NUCLEUS + 1},
  {0x0018, 0x0086, &Siemens_header.G18.Acq.EchoNumber, create_long_element, 1},
  {0x0018, 0x0090, &Siemens_header.G18.Acq.DataCollectionDiameter, create_long_element, 1},
  {0x0018, 0x1000, &Siemens_header.G18.Acq.DeviceSerialNumber, create_char_element, LENGTH_LABEL + 1},
  {0x0018, 0x1020, &Siemens_header.G18.Acq.SoftwareVersion, create_char_element, LENGTH_SOFTWARE_VERSION + 1},
  {0x0018, 0x1200, &Siemens_header.G18.Acq.CalibrationDate, create_ds_date_t_element, 1},
  {0x0018, 0x1201, &Siemens_header.G18.Acq.CalibrationTime, create_ds_time_t_element, 1},
  {0x0018, 0x1250, &Siemens_header.G18.Acq.ReceivingCoil, create_char_element, LENGTH_LABEL + 1},
  {0x0018, 0x5100, &Siemens_header.G18.Acq.PatientPosition, create_patient_position_t_element, 1},
  {0x0019, 0x1010, &Siemens_header.G19.Acq1.CM.NetFrequency, create_long_element, 1},
  {0x0019, 0x1020, &Siemens_header.G19.Acq1.CM.MeasurementMode, create_measurement_mode_t_element, 1},
  {0x0019, 0x1030, &Siemens_header.G19.Acq1.CM.CalculationMode, create_calculation_mode_t_element, 1},
  {0x0019, 0x1050, &Siemens_header.G19.Acq1.CM.NoiseLevel, create_long_element, 1},
  {0x0019, 0x1060, &Siemens_header.G19.Acq1.CM.NumberOfDataBytes, create_long_element, 1},
  {0x0019, 0x1210, &Siemens_header.G19.Acq2.Mr.TotalMeasurementTime, create_double_element, 1},
  {0x0019, 0x1211, &Siemens_header.G19.Acq2.Mr.TotalMeasurementTimeCur, create_double_element, 1},
  {0x0019, 0x1212, &Siemens_header.G19.Acq2.Mr.StartDelayTime, create_double_element, 1},
  {0x0019, 0x1213, &Siemens_header.G19.Acq2.Mr.DwellTime, create_double_element, 1},
  {0x0019, 0x1214, &Siemens_header.G19.Acq2.Mr.NumberOfPhases, create_long_element, 1},
  {0x0019, 0x1220, &Siemens_header.G19.Acq2.Mr.NumberOfFourierLinesNominal, create_long_element, 1},
  {0x0019, 0x1221, &Siemens_header.G19.Acq2.Mr.NumberOfFourierLinesCurrent, create_long_element, 1},
  {0x0019, 0x1226, &Siemens_header.G19.Acq2.Mr.NumberOfFourierLinesAfterZero, create_long_element, 1},
  {0x0019, 0x1228, &Siemens_header.G19.Acq2.Mr.FirstMeasuredFourierLine, create_long_element, 1},
  {0x0019, 0x1230, &Siemens_header.G19.Acq2.Mr.AcquisitionColumns, create_long_element, 1},
  {0x0019, 0x1231, &Siemens_header.G19.Acq2.Mr.ReconstructionColumns, create_long_element, 1},
  {0x0019, 0x1250, &Siemens_header.G19.Acq2.Mr.NumberOfAverages, create_long_element, 1},
  {0x0019, 0x1260, &Siemens_header.G19.Acq2.Mr.FlipAngle, create_double_element, 1},
  {0x0019, 0x1270, &Siemens_header.G19.Acq2.Mr.NumberOfPrescans, create_long_element, 1},
  {0x0019, 0x1281, &Siemens_header.G19.Acq2.Mr.FilterTypeRawData, create_filter_type_t_element, 1},
  {0x0019, 0x1282, &Siemens_header.G19.Acq2.Mr.FilterParameterRawData, create_filter_parameter_t_element, 1},
  {0x0019, 0x1283, &Siemens_header.G19.Acq2.Mr.FilterTypeImageData, create_filter_type_image_t_element, 1},
  {0x0019, 0x1285, &Siemens_header.G19.Acq2.Mr.FilterTypePhaseCorrection, create_filter_type_t_element, 1},
  {0x0019, 0x1290, &Siemens_header.G19.Acq2.Mr.NumberOfSaturationRegions, create_long_element, 1},
  {0x0019, 0x1294, &Siemens_header.G19.Acq2.Mr.ImageRotationAngle, create_double_element, 1},
  {0x0019, 0x1298, &Siemens_header.G19.Acq2.Mr.CoilPosition, create_image_location_t_element, 1},
  {0x0019, 0x1412, &Siemens_header.G19.Acq3.Mr.MagneticFieldStrength, create_double_element, 1},
  {0x0019, 0x1414, &Siemens_header.G19.Acq3.Mr.ADCVoltage, create_double_element, 1},
  {0x0019, 0x1416, &Siemens_header.G19.Acq3.Mr.ADCOffset, create_double_element, 2},
  {0x0019, 0x1420, &Siemens_header.G19.Acq3.Mr.TransmitterAmplitude, create_double_element, 1},
  {0x0019, 0x1421, &Siemens_header.G19.Acq3.Mr.NumberOfTransmitterAmplitudes, create_long_element, 1},
  {0x0019, 0x1422, &Siemens_header.G19.Acq3.Mr.TransmitterAttenuator, create_double_element, 1},
  {0x0019, 0x1424, &Siemens_header.G19.Acq3.Mr.TransmitterCalibration, create_double_element, 1},
  {0x0019, 0x1426, &Siemens_header.G19.Acq3.Mr.TransmitterReference, create_double_element, 1},
  {0x0019, 0x1450, &Siemens_header.G19.Acq3.Mr.ReceiverTotalGain, create_double_element, 1},
  {0x0019, 0x1451, &Siemens_header.G19.Acq3.Mr.ReceiverAmplifierGain, create_double_element, 1},
  {0x0019, 0x1452, &Siemens_header.G19.Acq3.Mr.ReceiverPreamplifierGain, create_double_element, 1},
  {0x0019, 0x1454, &Siemens_header.G19.Acq3.Mr.ReceiverCableAttenuation, create_double_element, 1},
  {0x0019, 0x1455, &Siemens_header.G19.Acq3.Mr.ReceiverReferenceGain, create_double_element, 1},
  {0x0019, 0x1456, &Siemens_header.G19.Acq3.Mr.ReceiverFilterFrequency, create_long_element, 1},
  {0x0019, 0x1460, &Siemens_header.G19.Acq3.Mr.ReconstructionScaleFactor, create_double_element, 1},
  {0x0019, 0x1462, &Siemens_header.G19.Acq3.Mr.ReferenceScaleFactor, create_double_element, 1},
  {0x0019, 0x1470, &Siemens_header.G19.Acq3.Mr.PhaseGradientAmplitude, create_double_element, 1},
  {0x0019, 0x1471, &Siemens_header.G19.Acq3.Mr.ReadoutGradientAmplitude, create_double_element, 1},
  {0x0019, 0x1472, &Siemens_header.G19.Acq3.Mr.SelectionGradientAmplitude, create_double_element, 1},
  {0x0019, 0x1480, &Siemens_header.G19.Acq3.Mr.GradientDelayTime, create_gradient_delay_time_t_element, 1},
  {0x0019, 0x1482, &Siemens_header.G19.Acq3.Mr.TotalGradientDelayTime, create_double_element, 1},
  {0x0019, 0x1490, &Siemens_header.G19.Acq3.Mr.SensitivityCorrectionLabel, create_char_element, LENGTH_LABEL + 1},
  {0x0019, 0x14a0, &Siemens_header.G19.Acq3.Mr.RfWatchdogMask, create_long_element, 1},
  {0x0019, 0x14a2, &Siemens_header.G19.Acq3.Mr.RfPowerErrorIndicator, create_double_element, 1},
  {0x0019, 0x14a5, &Siemens_header.G19.Acq3.Mr.SarWholeBody, create_sar_sed_t_element, 1},
  {0x0019, 0x14a6, &Siemens_header.G19.Acq3.Mr.Sed, create_sar_sed_t_element, 1},
  {0x0019, 0x14b0, &Siemens_header.G19.Acq3.Mr.AdjustmentStatusMask, create_long_element, 1},
  {0x0019, 0x1510, &Siemens_header.G19.Acq4.CM.ParameterFileName, create_char_element, LENGTH_FILE_NAME + 1},
  {0x0019, 0x1511, &Siemens_header.G19.Acq4.CM.SequenceFileName, create_char_element, LENGTH_FILE_NAME + 1},
  {0x0019, 0x1512, &Siemens_header.G19.Acq4.CM.SequenceFileOwner, create_char_element, LENGTH_SEQUENCE_INFO + 1},
  {0x0019, 0x1513, &Siemens_header.G19.Acq4.CM.SequenceDescription, create_char_element, LENGTH_SEQUENCE_INFO + 1},
  {0x0020, 0x0010, &Siemens_header.G20.Rel.Study, create_long_element, 1},
  {0x0020, 0x0012, &Siemens_header.G20.Rel.Acquisition, create_long_element, 1},
  {0x0020, 0x0013, &Siemens_header.G20.Rel.Image, create_long_element, 1},
  {0x0020, 0x0050, &Siemens_header.G20.Rel.Location, create_long_element, 1},
  {0x0020, 0x0070, &Siemens_header.G20.Rel.ImageGeometryType, create_geometry_t_element, 1},
  {0x0020, 0x1001, &Siemens_header.G20.Rel.AcquisitionsInSeries, create_long_element, 1},
  {0x0020, 0x1020, &Siemens_header.G20.Rel.Reference, create_reference_t_element, 1},
  {0x0021, 0x1011, &Siemens_header.G21.Rel1.CM.Target, create_target_point_t_element, 1},
  {0x0021, 0x1020, &Siemens_header.G21.Rel1.CM.RoiMask, create_short_element, 1},
  {0x0021, 0x1120, &Siemens_header.G21.Rel1.CM.FoV, create_field_of_view_t_element, 1},
  {0x0021, 0x1122, &Siemens_header.G21.Rel1.CM.ImageMagnificationFactor, create_double_element, 1},
  {0x0021, 0x1130, &Siemens_header.G21.Rel1.CM.ViewDirection, create_view_direction_t_element, 1},
  {0x0021, 0x1132, &Siemens_header.G21.Rel1.CM.RestDirection, create_rest_direction_t_element, 1},
  {0x0021, 0x1160, &Siemens_header.G21.Rel1.CM.ImagePosition, create_image_location_t_element, 1},
  {0x0021, 0x1161, &Siemens_header.G21.Rel1.CM.ImageNormal, create_image_location_t_element, 1},
  {0x0021, 0x1163, &Siemens_header.G21.Rel1.CM.ImageDistance, create_double_element, 1},
  {0x0021, 0x1165, &Siemens_header.G21.Rel1.CM.ImagePositioningHistoryMask, create_short_element, 1},
  {0x0021, 0x116a, &Siemens_header.G21.Rel1.CM.ImageRow, create_image_location_t_element, 1},
  {0x0021, 0x116b, &Siemens_header.G21.Rel1.CM.ImageColumn, create_image_location_t_element, 1},
  {0x0021, 0x1170, &Siemens_header.G21.Rel1.CM.PatientOrientationSet1, create_patient_orientation_t_element, 1},
  {0x0021, 0x1171, &Siemens_header.G21.Rel1.CM.PatientOrientationSet2, create_patient_orientation_t_element, 1},
  {0x0021, 0x1180, &Siemens_header.G21.Rel1.CM.StudyName, create_char_element, LENGTH_LABEL + 1},
  {0x0021, 0x1182, &Siemens_header.G21.Rel1.CM.StudyType, create_study_type_t_element, 1},
  {0x0021, 0x1322, &Siemens_header.G21.Rel2.Mr.PhaseCorRowRec, create_long_element, 1},
  {0x0021, 0x1324, &Siemens_header.G21.Rel2.Mr.PhaseCorColRec, create_long_element, 1},
  {0x0021, 0x1330, &Siemens_header.G21.Rel2.Mr.NumberOf3DRawPartNom, create_long_element, 1},
  {0x0021, 0x1331, &Siemens_header.G21.Rel2.Mr.NumberOf3DRawPartCur, create_long_element, 1},
  {0x0021, 0x1334, &Siemens_header.G21.Rel2.Mr.NumberOf3DImaPart, create_long_element, 1},
  {0x0021, 0x1336, &Siemens_header.G21.Rel2.Mr.Actual3DImaPartNumber, create_long_element, 1},
  {0x0021, 0x1339, &Siemens_header.G21.Rel2.Mr.SlabThickness, create_double_element, 1},
  {0x0021, 0x1340, &Siemens_header.G21.Rel2.Mr.NumberOfSlicesNom, create_long_element, 1},
  {0x0021, 0x1341, &Siemens_header.G21.Rel2.Mr.NumberOfSlicesCur, create_long_element, 1},
  {0x0021, 0x1342, &Siemens_header.G21.Rel2.Mr.CurrentSliceNumber, create_long_element, 1},
  {0x0021, 0x1343, &Siemens_header.G21.Rel2.Mr.CurrentGroupNumber, create_long_element, 1},
  {0x0021, 0x1344, &Siemens_header.G21.Rel2.Mr.CurrentSliceDistanceFactor, create_double_element, 1},
  {0x0021, 0x134f, &Siemens_header.G21.Rel2.Mr.OrderOfSlices, create_order_of_slices_t_element, 1},
  {0x0021, 0x1356, &Siemens_header.G21.Rel2.Mr.RepetitionTime, create_double_element, 1},
  {0x0021, 0x1370, &Siemens_header.G21.Rel2.Mr.NumberOfEchoes, create_long_element, 1},
  {0x0028, 0x0005, &Siemens_header.G28.Pre.ImageDimension, create_short_element, 1},
  {0x0028, 0x0010, &Siemens_header.G28.Pre.Rows, create_short_element, 1},
  {0x0028, 0x0011, &Siemens_header.G28.Pre.Columns, create_short_element, 1},
  {0x0028, 0x0030, &Siemens_header.G28.Pre.PixelSize, create_pixel_size_t_element, 1},
  {0x0028, 0x0040, &Siemens_header.G28.Pre.ImageFormat, create_image_format_t_element, 1},
  {0x0028, 0x0060, &Siemens_header.G28.Pre.CompressionCode, create_compression_code_t_element, 1},
  {0x0028, 0x0100, &Siemens_header.G28.Pre.BitsAllocated, create_short_element, 1},
  {0x0028, 0x0101, &Siemens_header.G28.Pre.BitsStored, create_short_element, 1},
  {0x0028, 0x0102, &Siemens_header.G28.Pre.HighBit, create_short_element, 1},
  {0x0028, 0x0103, &Siemens_header.G28.Pre.PixelRepresentation, create_short_element, 1},
  {0x0028, 0x1050, &Siemens_header.G28.Pre.WindowCenter, create_windows_t_element, 1},
  {0x0028, 0x1051, &Siemens_header.G28.Pre.WindowWidth, create_windows_t_element, 1},
  {0x0028, 0x1052, &Siemens_header.G28.Pre.RescaleIntercept, create_long_element, 1},
  {0x0028, 0x1053, &Siemens_header.G28.Pre.RescaleSlope, create_long_element, 1},
  {0x0029, 0x1110, &Siemens_header.G29.Pre.WindowStyle, create_window_style_t_element, 1},
  {0x0029, 0x1120, &Siemens_header.G29.Pre.PixelQualityCode, create_pixel_quality_code_t_element, 1},
  {0x0029, 0x1152, &Siemens_header.G29.Pre.SortCode, create_long_element, 1},
  {0, 0, NULL, NULL, 0}
};
#endif
