/**
 * @file  LayerSurface.h
 * @brief Layer data object for MRI volume.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2012/08/08 20:50:46 $
 *    $Revision: 1.54 $
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

#ifndef LayerSurface_h
#define LayerSurface_h

#include "LayerEditable.h"
#include "vtkSmartPointer.h"
#include <QList>

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
class LayerPropertySurface;
class FSSurface;
class LayerMRI;
class SurfaceOverlay;
class SurfaceAnnotation;
class SurfaceLabel;
class SurfaceROI;
class SurfaceSpline;

class LayerSurface : public LayerEditable
{
  Q_OBJECT
public:
  LayerSurface( LayerMRI* mri = NULL, QObject* parent = NULL );
  virtual ~LayerSurface();

  bool LoadSurfaceFromFile();
  bool LoadVectorFromFile();
  bool LoadCurvatureFromFile( const QString& filename );
  bool LoadOverlayFromFile( const QString& filename, const QString& fn_reg, bool bCorrelation );
  bool LoadGenericOverlayFromFile( const QString& filename, const QString& fn_reg );
  bool LoadCorrelationFromFile( const QString& filename );
  bool LoadAnnotationFromFile( const QString& filename );
  bool LoadLabelFromFile( const QString& filename );
  bool LoadSplineFromFile(const QString& filename);

  void Append2DProps( vtkRenderer* renderer, int nPlane );
  void Append3DProps( vtkRenderer* renderer, bool* bSliceVisibility = NULL );
  bool HasProp( vtkProp* prop );

  bool SaveSurface( const QString& filename );
  bool SaveSurface( );

  void SetSlicePositionToWorldCenter();

  inline LayerPropertySurface* GetProperty()
  {
    return (LayerPropertySurface*)mProperty;
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

  void SetPatchFileName( const QString& fn )
  {
    m_sPatchFilename = fn;
  }

  void SetVectorFileName( const QString& fn )
  {
    m_sVectorFilename = fn;
  }

  void SetTargetFileName( const QString& fn )
  {
    m_sTargetFilename = fn;
  }

  int GetActiveSurface();

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

  SurfaceOverlay* GetOverlay( const QString& name );

  int GetActiveOverlayIndex();

  SurfaceOverlay* GetActiveOverlay();

  void SetActiveOverlay( int nOverlay );

  void SetActiveOverlay( const QString& name );

  SurfaceOverlay* GetCorrelationOverlay();

  void CopyCorrelationOverlay(LayerSurface* surf);

  // annotation functions
  int GetNumberOfAnnotations();

  SurfaceAnnotation* GetAnnotation( int n );

  SurfaceAnnotation* GetAnnotation( const QString& name );

  int GetActiveAnnotationIndex();

  SurfaceAnnotation* GetActiveAnnotation();

  void SetActiveAnnotation( int n );

  void SetActiveAnnotation( const QString& name );

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

  SurfaceSpline* GetSpline()
  {
    return m_spline;
  }

  int GetHemisphere();

  void RepositionSurface( LayerMRI* mri, int nVertex, double value, int size, double sigma, int flags = 0 );
  void RepositionSurface( LayerMRI* mri, int nVertex, double* pos, int size, double sigma, int flags = 0 );
  void RepositionVertex( int nVertex, double* pos);
  bool SmoothSurface(int nMethod, int niters, double lambda, double k_cutoff);
  void RemoveIntersections();

  void Undo();
  bool HasUndo();

  void UpdateCorrelationOverlayAtVertex( int nVertex );
  void UpdateCorrelationOverlay();

  bool HasValidVolumeGeometry();

  int GetNumberOfVertices();

  SurfaceROI* GetSurfaceROI()
  {
    return m_roi;
  }

public slots:
  void SetActiveSurface( int nSurfaceType );
  void UpdateOverlay( bool bAskRedraw = true );
  void SetLoadAllSurfaces(bool bLoadAll)
  {
    m_bLoadAll = bLoadAll;
  }

  void SetActiveLabelColor(const QColor& c);
  void SetActiveLabelOutline(bool bOutline);

  void SetActiveAnnotationOutline(bool bOutline);

Q_SIGNALS:
  void SurfaceAnnotationAdded( SurfaceAnnotation* );
  void SurfaceLabelAdded( SurfaceLabel* );
  void SurfaceOverlayAdded( SurfaceOverlay* );
  void SurfaceOverlyDataUpdated();
  void SurfaceCurvatureLoaded();
  void SurfaceVectorLoaded();
  void ActiveSurfaceChanged( int n );
  void ActiveOverlayChanged( int n );
  void ActiveAnnotationChanged( int n );
  void ActiveLabelChanged( int n );

protected slots:
  void UpdateOpacity();
  void UpdateColorMap();
  void UpdateEdgeThickness();
  void UpdateVectorPointSize();
  void UpdateRenderMode();
  void UpdateVertexRender();
  void UpdateMeshRender();
  void UpdateActorPositions();
  void UpdateROIPosition(double dx, double dy, double dz);
  void UpdateVectorActor2D();

protected:
  void InitializeSurface();
  void InitializeActors();
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

  QString   m_sPatchFilename;
  QString   m_sVectorFilename;
  QString   m_sTargetFilename;

  vtkSmartPointer<vtkActor>   m_sliceActor2D[3];
  vtkSmartPointer<vtkActor>   m_sliceActor3D[3];
  vtkSmartPointer<vtkActor>   m_vectorActor2D[3];

  // vtkLODActor*  m_mainActor;
  vtkSmartPointer<vtkActor>   m_mainActor;
  vtkSmartPointer<vtkActor>   m_vectorActor;
  vtkSmartPointer<vtkActor>   m_vertexActor;
  vtkSmartPointer<vtkActor>   m_vertexActor2D[3];
  vtkSmartPointer<vtkActor>   m_wireframeActor;

  QList<SurfaceOverlay*>    m_overlays;
  int         m_nActiveOverlay;

  QList<SurfaceAnnotation*> m_annotations;
  int         m_nActiveAnnotation;

  QList<SurfaceLabel*>      m_labels;
  int         m_nActiveLabel;

  SurfaceROI*           m_roi;
  SurfaceSpline*        m_spline;

  bool        m_bUndoable;
  bool        m_bVector2DPendingUpdate;
  bool        m_bLoadAll;
};

#endif


