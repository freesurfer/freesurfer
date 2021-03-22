/**
 * @brief Base VolumeFilterDilate class.
 *
 */
/*
 * Original Author: Ruopeng Wang
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
 *
 */

#ifndef VolumeFilterDilate_h
#define VolumeFilterDilate_h

#include "VolumeFilter.h"

class LayerMRI;

class VolumeFilterDilate : public VolumeFilter
{
public:
  VolumeFilterDilate( LayerMRI* input = 0, LayerMRI* output = 0, QObject* parent = 0 );

  QString GetName()
  {
    return "Dilate";
  }

protected:
  bool Execute();
};

#endif


