#include <stdio.h>

#include "nii_dicom.h"
#include "dcm2niix_fswrapper.h"

struct TDCMopts dcm2niix_fswrapper::tdcmOpts;

/* These are the TDCMopts defaults set in dcm2niix
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
isIgnoreSeriesInstanceUID = false    // if true, d.seriesUidCrc = d.seriesNum;
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
filename = "%s_%p"    // seriesNum_protocol -f "%s_%p"
outdir = "..."
indir = "..."
pigzname = '\000'
optsname = "~/.dcm2nii.ini"
indirParent = "..."
imageComments = ""
seriesNumber = nnn 
numSeries = 0
 */

// set TDCMopts defaults, overwrite settings to output in mgz orientation
void dcm2niix_fswrapper::setOpts(const char* dcmindir, const char* niioutdir, bool createBIDS)
{
  memset(&tdcmOpts, 0, sizeof(tdcmOpts));
  setDefaultOpts(&tdcmOpts, NULL);

  if (dcmindir != NULL)
    strcpy(tdcmOpts.indir, dcmindir);
  if (niioutdir != NULL)
    strcpy(tdcmOpts.outdir, niioutdir);

  strcpy(tdcmOpts.filename, "%4s.%p");

  // set the options for freesurfer mgz orientation
  tdcmOpts.isRotate3DAcq = false;
  tdcmOpts.isFlipY = false;
  tdcmOpts.isIgnoreSeriesInstanceUID = true;
  tdcmOpts.isCreateBIDS = createBIDS;
  tdcmOpts.isGz = false;
  tdcmOpts.isForceStackSameSeries = 1; // merge 2D slice '-m y'
  tdcmOpts.isForceStackDCE = false;
  //tdcmOpts.isForceOnsetTimes = false;
}

// interface to isDICOMfile() in nii_dicom.cpp
bool dcm2niix_fswrapper::isDICOM(const char* file)
{
  return isDICOMfile(file);
}

/*
 * interface to nii_loadDirCore() to search all dicom files from the directory input file is in,
 * and convert dicom files with the same series as given file.
 */
int dcm2niix_fswrapper::dcm2NiiOneSeries(const char* dcmfile)
{
  // get seriesNo for given dicom file
  struct TDICOMdata tdicomData = readDICOM((char*)dcmfile);

  double seriesNo = (double)tdicomData.seriesUidCrc;
  if (tdcmOpts.isIgnoreSeriesInstanceUID)
    seriesNo = (double)tdicomData.seriesNum;

  // set TDCMopts to convert just one series
  tdcmOpts.seriesNumber[0] = seriesNo;
  tdcmOpts.numSeries = 1;

  return nii_loadDirCore(tdcmOpts.indir, &tdcmOpts);
}

// interface to nii_getMrifsStruct()
MRIFSSTRUCT* dcm2niix_fswrapper::getMrifsStruct(void)
{
  return nii_getMrifsStruct();
}

// return nifti header saved in MRIFSSTRUCT
nifti_1_header* dcm2niix_fswrapper::getNiiHeader(void)
{
  MRIFSSTRUCT* mrifsStruct = getMrifsStruct();
  return &mrifsStruct->hdr0;
}

// return image data saved in MRIFSSTRUCT
const unsigned char* dcm2niix_fswrapper::getMRIimg(void)
{
  MRIFSSTRUCT* mrifsStruct = getMrifsStruct();
  return mrifsStruct->imgM;
}

void dcm2niix_fswrapper::dicomDump(const char* dicomdir)
{
  strcpy(tdcmOpts.indir, dicomdir);
  tdcmOpts.isDumpNotConvert = true;
  nii_loadDirCore(tdcmOpts.indir, &tdcmOpts);

  return;

#if 0
  struct TSearchList nameList;
  int nConvertTotal = 0;
#if defined(_WIN64) || defined(_WIN32) || defined(USING_R)
  nameList.maxItems = 24000; // larger requires more memory, smaller more passes
#else //UNIX, not R
  nameList.maxItems = 96000; // larger requires more memory, smaller more passes
#endif  

  if ((is_fileNotDir(opts->indir)) && isExt(opts->indir, ".txt")) {
    nameList.str = (char **)malloc((nameList.maxItems + 1) * sizeof(char *)); //reserve one pointer (32 or 64 bits) per potential file
    nameList.numItems = 0;
    FILE *fp = fopen(opts->indir, "r"); //textDICOM
    if (fp == NULL)
      return EXIT_FAILURE;
    char dcmname[2048];
    while (fgets(dcmname, sizeof(dcmname), fp)) {
      int sz = (int)strlen(dcmname);
      if (sz > 0 && dcmname[sz - 1] == '\n')
	dcmname[sz - 1] = 0; //Unix LF
      if (sz > 1 && dcmname[sz - 2] == '\r')
	dcmname[sz - 2] = 0; //Windows CR/LF
      if ((!is_fileexists(dcmname)) || (!is_fileNotDir(dcmname))) { //<-this will accept meta data
	fclose(fp);
	printError("Problem with file '%s'\n", dcmname);
	return EXIT_FAILURE;
      }
      if (nameList.numItems < nameList.maxItems) {
	nameList.str[nameList.numItems] = (char *)malloc(strlen(dcmname) + 1);
	strcpy(nameList.str[nameList.numItems], dcmname);
      }
      nameList.numItems++;
    }
    fclose(fp);
    if (nameList.numItems >= nameList.maxItems) {
      printError("Too many file names in '%s'\n", opts->indir);
      return EXIT_FAILURE;
    }
    if (nameList.numItems < 1)
      return kEXIT_NO_VALID_FILES_FOUND;
    printMessage("Found %lu files in '%s'\n", nameList.numItems, opts->indir);
  } else {
    //1: find filenames of dicom files: up to two passes if we found more files than we allocated memory
    for (int i = 0; i < 2; i++) {
      nameList.str = (char **)malloc((nameList.maxItems + 1) * sizeof(char *)); //reserve one pointer (32 or 64 bits) per potential file
      nameList.numItems = 0;
      int ret = searchDirForDICOM(indir, &nameList, opts->dirSearchDepth, 0, opts);
      if (ret == EXIT_SUCCESS) //e.g. converted ECAT
	nConvertTotal++;
      if (nameList.numItems <= nameList.maxItems)
	break;
      freeNameList(nameList);
      nameList.maxItems = nameList.numItems + 1;
      //printMessage("Second pass required, found %ld images\n", nameList.numItems);
    }

    if (nameList.numItems < 1) {
      if ((opts->dirSearchDepth > 0) && (nConvertTotal < 1))
	printError("Unable to find any DICOM images in %s (or subfolders %d deep)\n", indir, opts->dirSearchDepth);
      else //keep silent for dirSearchDepth = 0 - presumably searching multiple folders
      {
      };
      free(nameList.str); //ignore compile warning - memory only freed on first of 2 passes
      if (nConvertTotal > 0) return EXIT_SUCCESS; //e.g. converted ECAT
	return kEXIT_NO_VALID_FILES_FOUND;
    }
  }
  
  size_t nDcm = nameList.numItems;
  printMessage("Found %lu DICOM file(s)\n", nameList.numItems); //includes images and other non-image DICOMs

  struct TDICOMdata *dcmList = (struct TDICOMdata *)malloc(nameList.numItems * sizeof(struct TDICOMdata));
  struct TDTI4D *dti4D = (struct TDTI4D *)malloc(sizeof(struct TDTI4D));
  struct TDCMprefs prefs;
  opts2Prefs(opts, &prefs);
  bool compressionWarning = false;
  bool convertError = false;
  bool isDcmExt = isExt(opts->filename, ".dcm"); // "%r.dcm" with multi-echo should generate "1.dcm", "1e2.dcm"
  if (isDcmExt)
    opts->filename[strlen(opts->filename) - 4] = 0; // "%s_%r.dcm" -> "%s_%r"
  //consider OpenMP
  // g++-9 -I. main_console.cpp nii_foreign.cpp nii_dicom.cpp jpg_0XC3.cpp ujpeg.cpp nifti1_io_core.cpp nii_ortho.cpp nii_dicom_batch.cpp -o dcm2niix -DmyDisableOpenJPEG -fopenmp

  for (int i = 0; i < (int)nDcm; i++) {
    if ((isExt(nameList.str[i], ".par")) && (isDICOMfile(nameList.str[i]) < 1)) {
      //strcpy(opts->indir, nameList.str[i]); //set to original file name, not path
      dcmList[i].converted2NII = 1;
      int ret = convert_parRec(nameList.str[i], *opts);
      if (ret == EXIT_SUCCESS)
	nConvertTotal++;
      else
	convertError = true;
      continue;
    }

    dcmList[i] = readDICOMx(nameList.str[i], &prefs, dti4D); //ignore compile warning - memory only freed on first of 2 passes
    //dcmList[i] = readDICOMv(nameList.str[i], opts->isVerbose, opts->compressFlag, dti4D); //ignore compile warning - memory only freed on first of 2 passes
    if (opts->isIgnoreSeriesInstanceUID)
      dcmList[i].seriesUidCrc = dcmList[i].seriesNum;
    //if (!dcmList[i].isValid) printf(">>>>Not a valid DICOM %s\n", nameList.str[i]);
    
    if ((dcmList[i].isValid) && ((dti4D->sliceOrder[0] >= 0) || (dcmList[i].CSA.numDti > 1))) { //4D dataset: dti4D arrays require huge amounts of RAM - write this immediately
      struct TDCMsort dcmSort[1];
      fillTDCMsort(dcmSort[0], i, dcmList[i]);
      //printMessage("***MGH_FREESURFER***: calling saveDcm2Nii() (%s:%s:%d)\n", __FILE__, __func__, __LINE__);
      dcmList[i].converted2NII = 1;
      int ret = saveDcm2Nii(1, dcmSort, dcmList, &nameList, *opts, dti4D);
      if (ret == EXIT_SUCCESS)
	nConvertTotal++;
      else
	convertError = true;
    }

    if ((dcmList[i].compressionScheme != kCompressNone) && (!compressionWarning) && (opts->compressFlag != kCompressNone)) {
      compressionWarning = true; //generate once per conversion rather than once per image
      printMessage("Image Decompression is new: please validate conversions\n");
    }
    if (opts->isProgress)
      progressPct = reportProgress(progressPct, kStage1Frac + (kStage2Frac * (float)i / (float)nDcm)); //proportion correct, 0..100
  }
#endif
}
