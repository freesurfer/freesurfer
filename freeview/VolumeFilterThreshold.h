/**
 * @brief Base VolumeFilterThreshold class.
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

#ifndef VolumeFilterThreshold_h
#define VolumeFilterThreshold_h

#include "VolumeFilter.h"

class LayerMRI;

class VolumeFilterThreshold : public VolumeFilter
{
public:
  VolumeFilterThreshold( LayerMRI* input = 0, LayerMRI* output = 0, QObject* parent = 0 );

  void SetReplaceIn( bool val )
  {
    m_bReplaceIn = val;
  }

  void SetReplaceOut( bool val )
  {
    m_bReplaceOut = val;
  }

  void SetThreshold( double* th )
  {
    m_dThreshold[0] = th[0];
    m_dThreshold[1] = th[1];
  }

  void SetInValue( double val)
  {
    m_dInValue = val;
  }

  void SetOutValue(double val)
  {
    m_dOutValue = val;
  }

  QString GetName()
  {
    return "Threshold";
  }

protected:
  bool Execute();

  bool    m_bReplaceIn;
  bool    m_bReplaceOut;
  double  m_dThreshold[2];
  double  m_dInValue;
  double  m_dOutValue;
};

#endif


