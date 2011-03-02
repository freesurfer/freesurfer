/**
 * @file  LayerMRI.h
 * @brief Layer class for MRI volume.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:02 $
 *    $Revision: 1.56 $
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
class vtkSimpleLabelEdgeFilter;
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
class vtkImageResample;
class vtkUnsignedCharArray;
class LayerPropertiesMRI;
class FSVolume;
class wxWindow;
class wxCommandEvent;
class BuildContourThread;
class Contour2D;
class SurfaceRegion;
class SurfaceRegionGroups;

class LayerMRI : public LayerVolumeBase
{
  friend class BuildContourThread;
  friend class VolumeCropper;
  friend class SurfaceRegionGroups;
  
public:
  LayerMRI( LayerMRI* ref );
  virtual ~LayerMRI();

//  bool LoadVolumeFromFile( std::string filename );
  bool LoadVolumeFromFile( wxWindow* wnd, wxCommandEvent& event );
  bool Create( LayerMRI* mri, bool bCopyVoxel, int data_type = -1 );

  virtual void Append2DProps( vtkRenderer* renderer, int nPlane );
  virtual void Append3DProps( vtkRenderer* renderer, bool* bPlaneVisibility = NULL );
  bool HasProp( vtkProp* prop );
  
  void Remove2DProps( vtkRenderer* render, int nPlane );

//  void SetSliceNumber( int* sliceNumber );
  void SetSlicePositionToWorldCenter();

  virtual double GetVoxelValue( double* pos );
  double GetVoxelValueByOriginalIndex( int i, int j, int k );
  
  virtual std::string GetLabelName( double value );

  void RASToOriginalIndex( const double* pos, int* n );
  void OriginalIndexToRAS( const int* n, double* pos );

  inline LayerPropertiesMRI* GetProperties()
  {
    return (LayerPropertiesMRI*)mProperties;
  }

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

  void NativeRASToTkReg( const double* pos_in, double* pos_out );
  
  void TkRegToNativeRAS( const double* pos_in, double* pos_out );
  
  int GetNumberOfFrames();

  void SetActiveFrame( int nFrame );

  void GetRASCenter( double* pt );
  
  bool IsTransformed();
  
  void SetReorient( bool bReorient );
  
  void SetSampleMethod( int nSampleMethod )
  {
    m_nSampleMethod = nSampleMethod;
  }

  void SetConform( bool bConform );
  
  bool GetVoxelValueRange( const double* pt0, const double* pt1, 
                           int nPlane, double* range_out );
  
  bool GetVoxelStatsRectangle( const double* pt0, const double* pt1, 
                           int nPlane, double* mean_out, double* sd_out = NULL, int* cnt_out = NULL );
  
  bool GetVoxelsOnLine( const double* pt0, const double* pt1, int nPlane, int*& indice_out, double*& value_out, int* cnt_out );
  
  void ResetWindowLevel();
  
  int GetDataType();
  
  COLOR_TABLE* GetEmbeddedColorTable();
  
  void SnapToVoxelCenter( const double* pt_in, double* pt_out );
  
  int GetBuildContourThreadID()
  {
    return m_nThreadID;
  }
  
  void RealizeContourActor();
  
  void GetCurrentLabelStats( int nPlane, float* label_out, int* count_out, float* area_out );
  
  vtkImageData* GetSliceImageData( int nPlane );
  
  bool FloodFillByContour2D( double* ras, Contour2D* c2d );
  
  virtual void SetModified();
  
  SurfaceRegion* CreateNewSurfaceRegion( double* pt );
  
  void AddSurfaceRegionLoopPoint( double* pt );
  
  void CloseSurfaceRegion();
  
  SurfaceRegion* SelectSurfaceRegion( double* pos );
  SurfaceRegion* SelectSurfaceRegion( int nId );
  
  SurfaceRegion* GetCurrentSurfaceRegion()
  {
    return m_currentSurfaceRegion;
  }
  
  bool DeleteCurrentSurfaceRegion();
  
  int GetNumberOfSurfaceRegions()
  {
    return m_surfaceRegions.size();
  }
  
  bool SaveAllSurfaceRegions( wxString& fn );
  
  bool LoadSurfaceRegions( wxString& fn );
  
  const char* GetOrientationString();
  
  void SetCroppingBounds( double* bounds );
  
  virtual void GetDisplayBounds( double* bounds );
  
  bool SaveRegistration( const char* filename );
  
  void GetLabelStats( LayerMRI* label, int nPlane,
                      std::vector<int>& id, 
                      std::vector<int>& number, 
                      std::vector<double>& mean, 
                      std::vector<double>& std );
  
  bool SaveContourToFile( const char* filename );
  
  SurfaceRegionGroups* GetSurfaceRegionGroups()
  {
    return m_surfaceRegionGroups;
  }
  
protected:
  virtual bool DoRotate( std::vector<RotationElement>& rotations, 
                       wxWindow* wnd, 
                       wxCommandEvent& event );
  virtual void DoTranslate( double* offset );
  virtual void DoScale( double* rscale, int nSampleMethod );
  virtual void DoRestore();
  
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
  
  void ResetSurfaceRegionIds();
  
  std::vector<int> GetVoxelIndicesBetweenPoints( int* n0, int* n1 );
  
  void UpdateLabelOutline();
  void UpdateUpSampleMethod();
  
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

  // Pipeline ------------------------------------------------------------
  vtkSmartPointer<vtkImageReslice>      mReslice[3];
  vtkSmartPointer<vtkImageMapToColors>  mColorMap[3];
  vtkSmartPointer<vtkSimpleLabelEdgeFilter>   mEdgeFilter[3];
  vtkSmartPointer<vtkImageResample>     mResample[3];

  FSVolume*   m_volumeSource;
  FSVolume*   m_volumeRef;
  bool    m_bResampleToRAS;
  bool    m_bReorient;
  int     m_nSampleMethod;
  bool    m_bConform;
  
  vtkImageActor*  m_sliceActor2D[3];
  vtkImageActor*  m_sliceActor3D[3];
  
  vtkActor*       m_glyphActor2D[3];
  vtkActor*       m_glyphActor3D[3];
  
  struct SegmentationActor
  {
    int id;
    vtkActor* actor;
  };
        
  std::vector<SegmentationActor>  m_segActors;              
  
  vtkSmartPointer<vtkActor>   m_actorContour;
  vtkSmartPointer<vtkVolume>  m_propVolume;
  
  int         m_nThreadID;
  vtkSmartPointer<vtkActor>       m_actorContourTemp;
  
  std::vector<SurfaceRegion*>     m_surfaceRegions;
  SurfaceRegion*                  m_currentSurfaceRegion;
  SurfaceRegionGroups*            m_surfaceRegionGroups;
  
  int         m_nOrientationIndex[3];
  
private:
  double**    private_buf1_3x3;
  double**    private_buf2_3x3;    
};

#endif


