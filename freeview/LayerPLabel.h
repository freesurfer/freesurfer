/**
 * @brief Layer class for P-Label volumes.
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


