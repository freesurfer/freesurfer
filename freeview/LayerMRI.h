/**
 * @file  LayerMRI.h
 * @brief Layer class for MRI volume.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2009/12/23 05:35:55 $
 *    $Revision: 1.29 $
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
#include <vector>

extern "C"
{
#include "colortab.h"
} 

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
class vtkPolyData;
class vtkPolyDataAlgorithm;
class vtkUnsignedCharArray;
class LayerPropertiesMRI;
class FSVolume;
class wxWindow;
class wxCommandEvent;
class BuildContourThread;

class LayerMRI : public LayerVolumeBase
{
  friend class BuildContourThread;
  
public:
  LayerMRI( LayerMRI* ref );
  virtual ~LayerMRI();

//  bool LoadVolumeFromFile( std::string filename );
  bool LoadVolumeFromFile( wxWindow* wnd, wxCommandEvent& event );
  bool Create( LayerMRI* mri, bool bCopyVoxel, int data_type = -1 );

  virtual void Append2DProps( vtkRenderer* renderer, int nPlane );
  virtual void Append3DProps( vtkRenderer* renderer, bool* bPlaneVisibility = NULL );
  bool HasProp( vtkProp* prop );

//  void SetSliceNumber( int* sliceNumber );
  void SetSlicePositionToWorldCenter();

  virtual double GetVoxelValue( double* pos );
  double GetVoxelValueByOriginalIndex( int i, int j, int k );
  
  virtual std::string GetLabelName( double value );

  void RASToOriginalIndex( const double* pos, int* n );
  void OriginalIndexToRAS( const int* n, double* pos );

  LayerPropertiesMRI* GetProperties();

  virtual void DoListenToMessage ( std::string const iMessage, void* iData, void* sender );

  virtual void SetVisible( bool bVisible = true );
  virtual bool IsVisible();
  
  virtual void UpdateVoxelValueRange( double dValue );

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

  void GetRASCenter( double* pt );
  
  virtual bool Rotate( std::vector<RotationElement>& rotations, 
                       wxWindow* wnd, 
                       wxCommandEvent& event );
  
  void SetReorient( bool bReorient );
  
  void SetSampleMethod( int nSampleMethod )
  {
    m_nSampleMethod = nSampleMethod;
  }

  bool GetVoxelValueRange( const double* pt0, const double* pt1, 
                           int nPlane, double* range_out );
  
  bool GetVoxelStats( const double* pt0, const double* pt1, 
                           int nPlane, double* mean_out, double* sd_out = NULL );
  
  void ResetWindowLevel();
  
  int GetDataType();
  
  COLOR_TABLE* GetEmbeddedColorTable();
  
  void SnagToVoxelCenter( const double* pt_in, double* pt_out );
  
  int GetBuildContourThreadID()
  {
    return m_nThreadID;
  }
  
  void RealizeContourActor();
  
protected:
  virtual void SetModified();

  void InitializeVolume();
  void InitializeActors();
  void UpdateOpacity();
  void UpdateResliceInterpolation();
  void UpdateTextureSmoothing();
  void UpdateContour( int nSegIndex = -1 );
  void UpdateContourActor( int nSegIndex );
  void UpdateContourColor();
  void ShowContour();
  void UpdateVolumeRendering();
  void UpdateVectorActor();
  void UpdateVectorActor( int nPlane, vtkImageData* imagedata );
  virtual void UpdateVectorActor( int nPlane );
  
  void UpdateTensorActor();
  void UpdateTensorActor( int nPlane, vtkImageData* imagedata = NULL );
  
  void BuildTensorGlyph( vtkImageData* imagedata,
                                 int i, int j, int k, 
                                 double* pt, double scale, 
                                 vtkPolyData* sourcepolydata,
                                 vtkUnsignedCharArray* scalars,
                                 vtkPolyDataAlgorithm* a);
  
  virtual void UpdateColorMap();

  virtual void OnSlicePositionChanged( int nPlane );

  LayerPropertiesMRI*      mProperties;
  // Pipeline ------------------------------------------------------------
  vtkSmartPointer<vtkImageReslice>   mReslice[3];
  vtkSmartPointer<vtkImageMapToColors>  mColorMap[3];

  FSVolume*   m_volumeSource;
  FSVolume*   m_volumeRef;
  bool    m_bResampleToRAS;
  bool    m_bReorient;
  int     m_nSampleMethod;

  vtkImageActor*  m_sliceActor2D[3];
  vtkImageActor*  m_sliceActor3D[3];
  
  vtkActor*       m_glyphActor2D[3];
  vtkActor*       m_glyphActor3D[3];
  
  struct SegmentationActor
  {
    int id;
    vtkActor* actor;
  };
        
  std::vector<SegmentationActor>   m_segActors;              
  
  vtkActor*   m_actorContour;
  vtkVolume*  m_propVolume;
  
  int         m_nThreadID;
  vtkSmartPointer<vtkActor> m_actorContourTemp;
  
private:
  double**    private_buf1_3x3;
  double**    private_buf2_3x3;    
};

#endif


