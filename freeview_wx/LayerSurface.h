/**
 * @file  LayerSurface.h
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

#ifndef LayerSurface_h
#define LayerSurface_h

#include "LayerEditable.h"
#include "vtkSmartPointer.h"
#include <string>
#include <vector>

class vtkImageReslice;
class vtkImageMapToColors;
class vtkTransform;
class vtkTexture;
class vtkPolyDataMapper;
class vtkActor;
class vtkLODActor;
class vtkImageActor;
class vtkImageData;
class vtkPlane;
class vtkDecimatePro;
class vtkProp;
class LayerPropertiesSurface;
class FSSurface;
class wxWindow;
class wxCommandEvent;
class LayerMRI;
class SurfaceOverlay;
class SurfaceAnnotation;
class SurfaceLabel;

class LayerSurface : public LayerEditable
{
public:
  LayerSurface( LayerMRI* mri = NULL );
  virtual ~LayerSurface();

  bool LoadSurfaceFromFile( wxWindow* wnd, wxCommandEvent& event );
  bool LoadVectorFromFile( wxWindow* wnd, wxCommandEvent& event );
  bool LoadCurvatureFromFile( const char* filename );
  bool LoadOverlayFromFile( const char* filename );
  bool LoadAnnotationFromFile( const char* filename );
  bool LoadLabelFromFile( const char* filename );

  void Append2DProps( vtkRenderer* renderer, int nPlane );
  void Append3DProps( vtkRenderer* renderer, bool* bSliceVisibility = NULL );
  bool HasProp( vtkProp* prop );
  
  bool SaveSurface( const char* filename, wxWindow* wnd, wxCommandEvent& event );
  bool SaveSurface( wxWindow* wnd, wxCommandEvent& event );

  void SetSlicePositionToWorldCenter();

  inline LayerPropertiesSurface* GetProperties()
  {
    return (LayerPropertiesSurface*)mProperties;
  }
  
  virtual void SetVisible( bool bVisible = true );
  virtual bool IsVisible();

  int GetVertexIndexAtRAS( double* ras, double* distance );

  int GetVertexIndexAtTarget( double* ras, double* distance );

  bool GetRASAtVertex       ( int nVertex, double* ras_out );
  bool GetSurfaceRASAtVertex( int nVertex, double* ras_out );

  bool GetTargetAtVertex( int nVertex, double* ras );
  
  void GetSurfaceRASAtTarget( double* pos_in, double* ras_out );
  
  void GetTargetAtSurfaceRAS( double* ras_in, double* pos_out );

  FSSurface* GetSourceSurface()
  {
    return m_surfaceSource;
  }

  void SetPatchFileName( const char* fn )
  {
    m_sPatchFilename = fn;
  }
  
  void SetVectorFileName( const char* fn )
  {
    m_sVectorFilename = fn;
  }

  void SetTargetFileName( const char* fn )
  {
    m_sTargetFilename = fn;
  }
  
  int GetActiveSurface();

  void SetActiveSurface( int nSurfaceType );

  int GetNumberOfVectorSets();

  int GetActiveVector();

  void SetActiveVector( int nVector );
  
  void GetVectorAtVertex( int nVertex, double* vec_out );
  
  void GetNormalAtVertex( int nVertex, double* vec_out );
  
  bool HasCurvature();
  
  void GetCurvatureRange( double* range );
  
  double GetCurvatureValue( int nVertex );
  
  // overlay functions
  bool HasOverlay();
  
  int GetNumberOfOverlays();
  
  SurfaceOverlay* GetOverlay( int n );
  
  SurfaceOverlay* GetOverlay( const char* name );
  
  int GetActiveOverlayIndex();
  
  SurfaceOverlay* GetActiveOverlay();
  
  void SetActiveOverlay( int nOverlay );
  
  void SetActiveOverlay( const char* name );
  
  void UpdateOverlay( bool bAskRedraw = false );
  
  // annotation functions
  int GetNumberOfAnnotations();
  
  SurfaceAnnotation* GetAnnotation( int n );
  
  SurfaceAnnotation* GetAnnotation( const char* name );
  
  int GetActiveAnnotationIndex();
  
  SurfaceAnnotation* GetActiveAnnotation();
  
  void SetActiveAnnotation( int n );
  
  void SetActiveAnnotation( const char* name );
  
  void UpdateAnnotation( bool bAskRedraw = false );
  
  // label functions
  int GetNumberOfLabels();
  
  SurfaceLabel* GetLabel( int n );
  
  SurfaceLabel* GetActiveLabel()
  {
    return ( m_nActiveLabel >= 0 ? m_labels[m_nActiveLabel] : NULL );
  }
  
  int GetActiveLabelIndex()
  {
    return m_nActiveLabel;
  }
  
  void SetActiveLabel( int n );
  
  LayerMRI* GetRefVolume()
  {
    return m_volumeRef;
  }
  
  int GetHemisphere();
  
  void RepositionSurface( LayerMRI* mri, int nVertex, double value, int size, double sigma );
  void RepositionSurface( LayerMRI* mri, int nVertex, double* pos, int size, double sigma );
  
  void Undo();
  bool HasUndo();
  
  bool LoadCorrelationFromFile( const char* filename );
  
  void UpdateCorrelationOverlayAtVertex( int nVertex );
  void UpdateCorrelationOverlay();
  
  bool HasVolumeGeometry();

protected:
  virtual void DoListenToMessage ( std::string const iMessage, void* iData, void* sender );

  void InitializeSurface();
  void InitializeActors();
  void UpdateOpacity();
  void UpdateColorMap();
  void UpdateEdgeThickness();
  void UpdateVectorPointSize();
  void UpdateRenderMode();
  void UpdateVertexRender();
  void UpdateMeshRender();
  void UpdateActorPositions();
  void UpdateVectorActor2D();
  void MapLabels( unsigned char* data, int nVertexCount );

  virtual void OnSlicePositionChanged( int nPlane );

  // Pipeline ------------------------------------------------------------
  vtkSmartPointer<vtkPlane>     mReslicePlane[3];
  vtkSmartPointer<vtkImageMapToColors>  mColorMap[3];

  vtkSmartPointer<vtkDecimatePro>   mLowResFilter;
  vtkSmartPointer<vtkDecimatePro>   mMediumResFilter;

  FSSurface*   m_surfaceSource;
  bool    m_bResampleToRAS;
  LayerMRI*   m_volumeRef;

  std::string   m_sPatchFilename;
  std::string   m_sVectorFilename;
  std::string   m_sTargetFilename;

  vtkActor*   m_sliceActor2D[3];
  vtkActor*   m_sliceActor3D[3];
  vtkActor*   m_vectorActor2D[3];

  // vtkLODActor*  m_mainActor;
  vtkSmartPointer<vtkActor>   m_mainActor;
  vtkSmartPointer<vtkActor>   m_vectorActor;
  vtkSmartPointer<vtkActor>   m_vertexActor;
  vtkSmartPointer<vtkActor>   m_wireframeActor;
  
  std::vector<SurfaceOverlay*>    m_overlays;
  int         m_nActiveOverlay;
  
  std::vector<SurfaceAnnotation*> m_annotations;
  int         m_nActiveAnnotation;
  
  std::vector<SurfaceLabel*>      m_labels;
  int         m_nActiveLabel;
  
  bool        m_bUndoable;
  bool        m_bVector2DPendingUpdate;
};

#endif


