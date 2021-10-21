#ifndef DCM2FSWRAPPER_H
#define DCM2FSWRAPPER_H

#include "nii_dicom_batch.h"
#include "nii_dicom.h"

class dcm2fsWrapper
{
public:
  static void setOpts(const char* dcmindir, const char* niioutdir);
  static bool isDICOM(const char* file);
  static int dcm2NiiOneSeries(const char* dcmfile);
  static int saveNii(const char* nii);

  static MRIFSSTRUCT* getMrifsStruct(void);
  static nifti_1_header* getNiiHeader(void);
  static const unsigned char* getMRIimg(void);

private:
  static struct TDCMopts tdcmOpts;

  //static MRIFSSTRUCT mrifsStruct;
};

#endif
