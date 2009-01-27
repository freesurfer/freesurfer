/**
 * @file  LayerMRI.h
 * @brief Layer data object for MRI volume.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2009/01/27 18:43:48 $
 *    $Revision: 1.3.2.2 $
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

#ifndef LayerMRI_h
#define LayerMRI_h

#include "LayerVolumeBase.h"
#include "vtkSmartPointer.h"
#include "CommonDataStruct.h"
#include <string>

class vtkFSVolumeSource;
class vtkImageReslice;
class vtkImageMapToColors;
class vtkTransform;
class vtkTexture;
class vtkPolyDataMapper;
class vtkActor;
class vtkImageActor;
class vtkImageData;
class vtkProp;
class vtkVolume;
class LayerPropertiesMRI;
class FSVolume;
class wxWindow;
class wxCommandEvent;

class LayerMRI : public LayerVolumeBase
{
public:
  LayerMRI( LayerMRI* ref );
  virtual ~LayerMRI();

//  bool LoadVolumeFromFile( std::string filename );
  bool LoadVolumeFromFile( wxWindow* wnd, wxCommandEvent& event );
  bool Create( LayerMRI* mri, bool bCopyVoxel );

  void Append2DProps( vtkRenderer* renderer, int nPlane );
  void Append3DProps( vtkRenderer* renderer, bool* bPlaneVisibility = NULL );
  bool HasProp( vtkProp* prop );

//  void SetSliceNumber( int* sliceNumber );
  void SetSlicePositionToWorldCenter();

  virtual double GetVoxelValue( double* pos );
  virtual std::string GetLabelName( double value );

  void RASToOriginalIndex( const double* pos, int* n );
  void OriginalIndexToRAS( const int* n, double* pos );

  LayerPropertiesMRI* GetProperties();

  virtual void DoListenToMessage ( std::string const iMessage, void* const iData );

  virtual void SetVisible( bool bVisible = true );
  virtual bool IsVisible();

  FSVolume* GetSourceVolume()
  {
    return m_volumeSource;
  }

  FSVolume* GetRefVolume()
  {
    return m_volumeRef;
  }

  void SetRefVolume( FSVolume* ref )
  {
    m_volumeRef = ref;
  }

  bool SaveVolume( wxWindow* wnd, wxCommandEvent& event );

  void SetResampleToRAS( bool bResample );

  bool GetResampleToRAS()
  {
    return m_bResampleToRAS;
  }

  void RemapPositionToRealRAS( const double* pos_in, double* pos_out );

  void TargetToRAS( const double* pos_in, double* pos_out )
  {
    RemapPositionToRealRAS( pos_in, pos_out );
  }

  void RemapPositionToRealRAS( double x_in, double y_in, double z_in,
                               double& x_out, double& y_out, double& z_out );
  void RASToTarget( const double* pos_in, double* pos_out );

  int GetNumberOfFrames();

  void SetActiveFrame( int nFrame );

  virtual bool Rotate( std::vector<RotationElement>& rotations, wxWindow* wnd, wxCommandEvent& event );

protected:
  virtual void SetModified();

  void InitializeVolume();
  void InitializeActors();
  void UpdateOpacity();
  void UpdateResliceInterpolation();
  void UpdateTextureSmoothing();
  void UpdateContour();
  void UpdateVolumeRendering();
  virtual void UpdateColorMap();

  virtual void OnSlicePositionChanged( int nPlane );

  LayerPropertiesMRI*      mProperties;
  // Pipeline ------------------------------------------------------------
  vtkSmartPointer<vtkImageReslice>   mReslice[3];
  vtkSmartPointer<vtkImageMapToColors>  mColorMap[3];

  FSVolume*   m_volumeSource;
  FSVolume*   m_volumeRef;
  bool    m_bResampleToRAS;

  vtkImageActor*  m_sliceActor2D[3];
  vtkImageActor*  m_sliceActor3D[3];
  vtkActor*   m_actorContour;
  vtkVolume*   m_propVolume;
};

#endif


