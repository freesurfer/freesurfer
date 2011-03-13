/**
 * @file  VolumeFilterMedian.h
 * @brief Base VolumeFilterMedian class.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/13 23:04:18 $
 *    $Revision: 1.4 $
 *
 * Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
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

#ifndef VolumeFilterMedian_h
#define VolumeFilterMedian_h

#include "VolumeFilter.h"

class LayerMRI;

class VolumeFilterMedian : public VolumeFilter
{
public:
  VolumeFilterMedian( LayerMRI* input = 0, LayerMRI* output = 0, QObject* parent = 0 );

  QString GetName()
  {
    return "Median";
  }

protected:
  bool Execute();
};

#endif


