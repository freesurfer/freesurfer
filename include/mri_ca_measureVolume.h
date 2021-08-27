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


#ifndef MeasureVolume_h
#define MeasureVolume_h

#include <string>


#include "mri.h"


using namespace std;




class CMeasureVolume
{
public:
  enum MeasureType { T1, T2, PD, EEG, MEG, fMRI, Unknown};


  MeasureType measureType;
  MRI* pVolume;
  string strMeasureFileDir;
};


#endif
