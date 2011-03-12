/**
 * @file  LayerPLabel.h
 * @brief Layer class for P-Label volumes.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: krish $
 *    $Date: 2011/03/12 00:28:49 $
 *    $Revision: 1.5 $
 *
 * Copyright (C) 2008-2009,
 * The General Hospital Corporation (Boston, MA).
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 * Bug reports: analysis-bugs@nmr.mgh.harvard.edu
 *
 */

#ifndef LayerPLabel_h
#define LayerPLabel_h

#include "LayerMRI.h"
#include "vtkSmartPointer.h"
#include <QStringList>

class FSVolume;
class vtkImageData;

class LayerPLabel : public LayerMRI
{
    Q_OBJECT
public:
  LayerPLabel( LayerMRI* ref, QObject* parent = NULL );
  virtual ~LayerPLabel();

  bool LoadVolumeFiles();

  void SetVolumeFileNames( const QStringList& filenames )
  {
    m_sFilenames = filenames;
  }

  void SetFileNamePrefix( const QString& prefix )
  {
    m_sFilenamePrefix = prefix;
  }

  void SetLUT( const QString& lut )
  {
    m_sLUT = lut;
  }

  double GetVoxelValue(double* pos);

  QString GetLabelName(double* pos);

protected:
//  bool DoRotate( std::vector<RotationElement>& rotations, wxWindow* wnd, wxCommandEvent& event );
  void UpdateColorMap();

  FSVolume*     m_volumeTemp;
  QStringList   m_sFilenames;
  QString       m_sFilenamePrefix;
  QString       m_sLUT;
  vtkSmartPointer<vtkImageData> m_imageIndex;
};

#endif


