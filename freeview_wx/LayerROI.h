/**
 * @file  LayerROI.h
 * @brief Layer data object for MRI volume.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: krish $
 *    $Date: 2011/03/11 23:27:39 $
 *    $Revision: 1.1 $
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
 */

#ifndef LayerROI_h
#define LayerROI_h

#include "LayerVolumeBase.h"
#include "vtkSmartPointer.h"
#include <string>
#include <vector>

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
class LayerPropertiesROI;
class wxWindow;
class wxCommandEvent;

class LayerROI : public LayerVolumeBase
{
public:
  LayerROI( LayerMRI* layerMRI );
  virtual ~LayerROI();

  bool LoadROIFromFile( std::string filename );

  virtual void Append2DProps( vtkRenderer* renderer, int nPlane );
  virtual void Append3DProps( vtkRenderer* renderer, bool* bSliceVisibility = NULL );

  bool HasProp( vtkProp* prop );

  virtual void DoListenToMessage ( std::string const iMessage, void* iData, void* sender );

  void SetVisible( bool bVisible = true );
  bool IsVisible();

  inline LayerPropertiesROI* GetProperties()
  {
    return (LayerPropertiesROI*)mProperties;
  }

  bool SaveROI( wxWindow* wnd, wxCommandEvent& event );

  void UpdateLabelData( wxWindow* wnd, wxCommandEvent& event );

  virtual void SetModified();
  
protected:
  bool DoRotate( std::vector<RotationElement>& rotations, wxWindow* wnd, wxCommandEvent& event );
  void DoRestore();
  void InitializeActors();
  void UpdateOpacity();
  void UpdateColorMap();

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


