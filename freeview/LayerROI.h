/**
 * @file  LayerROI.h
 * @brief Layer data object for MRI volume.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2009/01/27 18:43:48 $
 *    $Revision: 1.2.2.2 $
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

  virtual void DoListenToMessage ( std::string const iMessage, void* const iData );

  void SetVisible( bool bVisible = true );
  bool IsVisible();

  LayerPropertiesROI* GetProperties();

  bool SaveROI( wxWindow* wnd, wxCommandEvent& event );

  void UpdateLabelData( wxWindow* wnd, wxCommandEvent& event );

  bool Rotate( std::vector<RotationElement>& rotations, wxWindow* wnd, wxCommandEvent& event );

protected:
  void InitializeActors();
  void UpdateOpacity();
  void UpdateColorMap();
  virtual void SetModified();

  virtual void OnSlicePositionChanged( int nPlane );

  LayerPropertiesROI*      mROIProperties;

  // Pipeline ------------------------------------------------------------
  vtkSmartPointer<vtkImageReslice>   mReslice[3];
  vtkSmartPointer<vtkImageMapToColors>  mColorMap[3];

  LayerMRI*   m_layerSource;
  FSLabel*   m_label;

  vtkImageActor*  m_sliceActor2D[3];
  vtkImageActor*  m_sliceActor3D[3];
};

#endif


