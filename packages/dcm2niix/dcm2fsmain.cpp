#include <stdio.h>
#include <stdlib.h>

#include "dcm2fsWrapper.h"

int main(int argc, const char *argv[])
{
  const char* dcmfile = argv[1];
  const char* niifile = argv[2];

#if 1  
  if (!dcm2fsWrapper::isDICOM(dcmfile)) {
    printf("ERROR: %s is not a dicom file or some other problem\n", dcmfile);
    exit(1);
  }

  dcm2fsWrapper::setOpts(NULL, NULL);
  dcm2fsWrapper::dcm2NiiOneSeries(dcmfile);
  dcm2fsWrapper::saveNii(niifile);

  return 0;
#else
  struct TDICOMdata tdicomData = readDICOM((char*)dcmfile);
  int seriesNo = tdicomData.seriesNum;

  // set the TDCMopts
  struct TDCMopts tdcmOpts;
  setDefaultOpts(&tdcmOpts, NULL);
  tdcmOpts.seriesNumber[0] = seriesNo;
  tdcmOpts.numSeries = 1;
  tdcmOpts.isIgnoreSeriesInstanceUID = true;
  tdcmOpts.isCreateBIDS = false;
  strcpy(tdcmOpts.indir, "/space/papancha/1/users/yh887/testdata/tutorial_data_20190918_1558/practice_with_data/DICOM");
  strcpy(tdcmOpts.outdir, "/space/papancha/1/users/yh887/dcm2niix_out");
  //getFileNameX(tdcmOpts.indir, dcmfile, 512);

 

  NIFTISTRUCT niftiStruct;
  nii_loadDirCore(tdcmOpts.indir, &tdcmOpts, &niftiStruct);

  FILE *fp = fopen("/space/papancha/1/users/yh887/dcm2niix_out/fs.nii", "wb");
  if (!fp) 
    return 1;

  //if (!opts.isSaveNativeEndian) swapEndian(&hdr, im, true); //byte-swap endian (e.g. little->big)
  fwrite(&niftiStruct.hdr0, sizeof(niftiStruct.hdr0), 1, fp);
  uint32_t pad = 0;
  fwrite(&pad, sizeof( pad), 1, fp);
  fwrite(niftiStruct.imgM, niftiStruct.imgsz, 1, fp);
  fclose(fp);

  return 0;
#endif
}
