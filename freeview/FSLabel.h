/**
 * @file  FSLabel.h
 * @brief Base label class that takes care of I/O and data conversion.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2016/12/06 18:25:54 $
 *    $Revision: 1.20 $
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
class FSSurface;

class FSLabel : public QObject
{
  Q_OBJECT
public:
  FSLabel( QObject* parent );
  virtual ~FSLabel();

  bool LabelRead( const QString& filename );
  bool LabelWrite( const QString& filename );

  void Initialize(FSVolume* ref_vol, FSSurface* surf, int coords);

  void UpdateLabelFromImage( vtkImageData* rasImage_in, FSVolume* ref_vol );
  void UpdateRASImage( vtkImageData* rasImage_out, FSVolume* ref_vol, double threshold = -1e10 );
  void FillUnassignedVertices(FSSurface* surf, FSVolume* mri_template, int coords);
  void EditVoxel(int nx, int ny, int nz, bool bAdd);

  bool GetCentroidRASPosition(double* pos, FSVolume* ref_vol);

  void GetStatsRange(double* range);

  LABEL* GetRawLabel()
  {
      return m_label;
  }

protected:
  LABEL*   m_label;
  double   m_dStatsRange[2];
};

#endif


