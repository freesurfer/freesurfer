/**
 * @brief Base VolumeFilterGradient class.
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

#ifndef VolumeFilterGradient_h
#define VolumeFilterGradient_h

#include "VolumeFilter.h"

class LayerMRI;

class VolumeFilterGradient : public VolumeFilter
{
public:
  VolumeFilterGradient( LayerMRI* input = 0, LayerMRI* output = 0, QObject* parent = 0 );

  void SetSmoothing( bool bSmooth )
  {
    m_bSmoothing = bSmooth;
  }

  bool GetSmoothing()
  {
    return m_bSmoothing;
  }

  void SetStandardDeviation( double sd )
  {
    m_dSD = sd;
  }

  double GetStandardDeviation()
  {
    return m_dSD;
  }

  QString GetName()
  {
    return "Gradient";
  }

protected:
  bool Execute();

  bool    m_bSmoothing;
  double  m_dSD;
};

#endif


