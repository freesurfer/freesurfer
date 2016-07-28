/**
 * @file  FSLabel.h
 * @brief Base label class that takes care of I/O and data conversion.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: zkaufman $
 *    $Date: 2016/07/28 14:52:37 $
 *    $Revision: 1.17.2.1 $
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

#ifndef FSLabel_h
#define FSLabel_h

#include <QObject>
#include "vtkImageData.h"
#include "vtkMatrix4x4.h"

extern "C"
{
#include "label.h"
}

class FSVolume;

class FSLabel : public QObject
{
  Q_OBJECT
public:
  FSLabel( QObject* parent );
  virtual ~FSLabel();

  bool LabelRead( const QString& filename );
  bool LabelWrite( const QString& filename );

  void UpdateLabelFromImage( vtkImageData* rasImage_in, FSVolume* ref_vol );
  void UpdateRASImage( vtkImageData* rasImage_out, FSVolume* ref_vol, double threshold = -1e10 );

  bool GetCentroidRASPosition(double* pos, FSVolume* ref_vol);

  void GetStatsRange(double* range);

protected:
  LABEL*   m_label;
  double   m_dStatsRange[2];
};

#endif


