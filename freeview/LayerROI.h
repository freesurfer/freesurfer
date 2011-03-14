/**
 * @file  LayerROI.h
 * @brief Layer data object for MRI volume.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/14 23:44:47 $
 *    $Revision: 1.17 $
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

#ifndef LayerROI_h
#define LayerROI_h

#include "LayerVolumeBase.h"
#include "vtkSmartPointer.h"

class FSLabel;
class vtkImageReslice;
class vtkImageMapToColors;
class vtkTransform;
class vtkTexture;
class vtkPolyDataMapper;
class vtkActor;
class vtkImageActor;
class vtkImageData;
class vtkProp;
class LayerMRI;
class LayerPropertyROI;

class LayerROI : public LayerVolumeBase
{
  Q_OBJECT
public:
  LayerROI( LayerMRI* layerMRI, QObject* parent = NULL  );
  virtual ~LayerROI();

  bool LoadROIFromFile( const QString& filename );

  virtual void Append2DProps( vtkRenderer* renderer, int nPlane );
  virtual void Append3DProps( vtkRenderer* renderer, bool* bSliceVisibility = NULL );

  bool HasProp( vtkProp* prop );

  void SetVisible( bool bVisible = true );
  bool IsVisible();

  inline LayerPropertyROI* GetProperty()
  {
    return (LayerPropertyROI*)mProperty;
  }

  bool SaveROI();

  void UpdateLabelData();

  virtual void SetModified();

protected slots:
  void UpdateOpacity();
  void UpdateColorMap();

protected:
  bool DoRotate( std::vector<RotationElement>& rotations );
  void DoRestore();
  void InitializeActors();

  virtual void OnSlicePositionChanged( int nPlane );

  // Pipeline ------------------------------------------------------------
  vtkSmartPointer<vtkImageReslice>   mReslice[3];
  vtkSmartPointer<vtkImageMapToColors>  mColorMap[3];

  LayerMRI*   m_layerSource;
  FSLabel*   m_label;

  vtkImageActor*  m_sliceActor2D[3];
  vtkImageActor*  m_sliceActor3D[3];
};

#endif


