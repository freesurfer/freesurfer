#include <stdio.h>

#include "nii_dicom.h"
#include "dcm2fsWrapper.h"

struct TDCMopts dcm2fsWrapper::tdcmOpts;
//MRIFSSTRUCT dcm2fsWrapper::mrifsStruct;

/* oct06 version
isIgnoreTriggerTimes = false
isTestx0021x105E = false
isAddNamePostFixes = true
isSaveNativeEndian = true
isOneDirAtATime = false
isRenameNotConvert = false
isSave3D = false
isGz = false
isPipedGz = false
isFlipY = true
isCreateBIDS = true
isSortDTIbyBVal = false
isAnonymizeBIDS = true
isOnlyBIDS = false
isCreateText = false
isForceOnsetTimes = true
isIgnoreDerivedAnd2D = false
isPhilipsFloatNotDisplayScaling = true
isTiltCorrect = true
isRGBplanar = false
isOnlySingleFile = false
isForceStackDCE = true
isIgnoreSeriesInstanceUID = false    // if true d.seriesUidCrc = d.seriesNum;
isRotate3DAcq = true
isCrop = false
saveFormat = 0
isMaximize16BitRange = 2
isForceStackSameSeries = 2
nameConflictBehavior = 2
isVerbose = 0
isProgress = 0
compressFlag = 0
dirSearchDepth = 5
gzLevel = 6
filename = "%s_%p\000%t_%s"
outdir = "/space/papancha/1/users/yh887/fs_test/dcm2niix/test2"
indir = "/autofs/space/sulc_001/users/xdicom/dicom/INTEROPERABILITY"
pigzname = '\000'
optsname = "/homes/7/yh887/.dcm2nii.ini"
indirParent = "INTEROPERABILITY", '\000' <repeats 495 times>, 
imageComments = "\000\325\377\377\377\177\000\000`\325\377\377\377\177\000\000&\260be\000\000\000", 
seriesNumber = {0 <repeats 16 times>}, 
numSeries = 0
 */
void dcm2fsWrapper::setOpts(const char* dcmindir, const char* niioutdir)
{
  memset(&tdcmOpts, 0, sizeof(tdcmOpts));
  setDefaultOpts(&tdcmOpts, NULL);

  if (dcmindir != NULL)
    strcpy(tdcmOpts.indir, dcmindir);
  if (niioutdir != NULL)
    strcpy(tdcmOpts.outdir, niioutdir);

  tdcmOpts.isRotate3DAcq = false;
  tdcmOpts.isFlipY = false;
  tdcmOpts.isIgnoreSeriesInstanceUID = true;
  tdcmOpts.isCreateBIDS = false;
  tdcmOpts.isGz = false;
  //tdcmOpts.isForceStackSameSeries = 1; // merge 2D slice '-m y'
  tdcmOpts.isForceStackDCE = false;
  //tdcmOpts.isForceOnsetTimes = false;
}

bool dcm2fsWrapper::isDICOM(const char* file)
{
  return isDICOMfile(file);
}

int dcm2fsWrapper::dcm2NiiOneSeries(const char* dcmfile)
{
  struct TDICOMdata tdicomData = readDICOM((char*)dcmfile);
  double seriesNo = (double)tdicomData.seriesUidCrc;  //seriesNum;
  tdcmOpts.seriesNumber[0] = seriesNo;
  tdcmOpts.numSeries = 1;

  //memset(&mrifsStruct, 0, sizeof(mrifsStruct));
  return nii_loadDirCore(tdcmOpts.indir, &tdcmOpts);
}

int dcm2fsWrapper::saveNii(const char* nii)
{
  FILE *fp = fopen(nii, "wb");
  if (!fp) 
    return 1;

  MRIFSSTRUCT* mrifsStruct = getMrifsStruct();

  //if (!opts.isSaveNativeEndian) swapEndian(&hdr, im, true); //byte-swap endian (e.g. little->big)
  fwrite(&mrifsStruct->hdr0, sizeof(mrifsStruct->hdr0), 1, fp);
  uint32_t pad = 0;
  fwrite(&pad, sizeof( pad), 1, fp);
  fwrite(mrifsStruct->imgM, mrifsStruct->imgsz, 1, fp);
  fclose(fp);

  return 0;
}

MRIFSSTRUCT* dcm2fsWrapper::getMrifsStruct(void)
{
  return nii_getMrifsStruct();
}

nifti_1_header* dcm2fsWrapper::getNiiHeader(void)
{
  MRIFSSTRUCT* mrifsStruct = getMrifsStruct();
  return &mrifsStruct->hdr0;
}

const unsigned char* dcm2fsWrapper::getMRIimg(void)
{
  MRIFSSTRUCT* mrifsStruct = getMrifsStruct();
  return mrifsStruct->imgM;
}
