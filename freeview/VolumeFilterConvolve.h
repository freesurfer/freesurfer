/**
 * @brief Base VolumeFilterConvolve class.
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

#ifndef VolumeFilterConvolve_h
#define VolumeFilterConvolve_h

#include "VolumeFilter.h"

class LayerMRI;

class VolumeFilterConvolve : public VolumeFilter
{
public:
  VolumeFilterConvolve( LayerMRI* input = 0, LayerMRI* output = 0, QObject* parent = 0 );

  void SetSigma( double dvalue )
  {
    m_dSigma = dvalue;
  }

  double GetSigma()
  {
    return m_dSigma;
  }

  QString GetName()
  {
    return "Convolve";
  }

protected:
  bool Execute();

  double m_dSigma;
};

#endif


