/**
 * @brief Layer data object for MRI volume.
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

#ifndef LayerSurface_h
#define LayerSurface_h

#include "LayerEditable.h"
#include "vtkSmartPointer.h"
#include <QList>
#include <QVector>
#include <QVariantMap>

class vtkImageReslice;
class vtkImageMapToColors;
class vtkTransform;
class vtkTexture;
class vtkPolyDataMapper;
class vtkActor;
class vtkImageActor;
class vtkImageData;
class vtkPlane;
class vtkDecimatePro;
class vtkProp;
class vtkCutter;
class vtkBox;
class LayerPropertySurface;
class FSSurface;
class LayerMRI;
class SurfaceOverlay;
class SurfaceAnnotation;
class SurfaceLabel;
class SurfaceROI;
class SurfaceSpline;
class SurfacePath;
class LayerROI;
class vtkPolyData;

struct RGBMap {
  QString name;
  QList<int> data;
};

class LayerSurface : public LayerEditable
{
  Q_OBJECT
public:
  LayerSurface( LayerMRI* mri = NULL, QObject* parent = NULL );
  virtual ~LayerSurface();

  bool LoadSurfaceFromFile(bool bIgnoreVG, QString& sAffineXformFilename);
  bool LoadVectorFromFile();
  bool LoadCurvatureFromFile( const QString& filename );
  bool LoadOverlayFromFile( const QString& filename, const QString& fn_reg, bool bCorrelation, bool bSecondHalfData = false );
  bool LoadGenericOverlayFromFile( const QString& filename, const QString& fn_reg, bool bSecondHalfData = false );
  bool LoadCorrelationFromFile( const QString& filename );
  bool LoadAnnotationFromFile( const QString& filename );
  bool LoadLabelFromFile( const QString& filename );
  bool LoadSplineFromFile(const QString& filename);
  bool LoadRGBFromFile(const QString& filename);
  bool CreateFromMRIS(void* mris_ptr);

  void Append2DProps( vtkRenderer* renderer, int nPlane );
  void Append3DProps( vtkRenderer* renderer, bool* bSliceVisibility = NULL );
  bool HasProp( vtkProp* prop );

  bool SaveSurface( const QString& filename );
  bool SaveSurface( );
  bool WriteIntersection( const QString& filename, int nPlane, LayerMRI* ref);
  bool SaveSurfaceAsSTL(const QString& filename);

  void SetSlicePositionToWorldCenter();

  inline LayerPropertySurface* GetProperty()
  {
    return (LayerPropertySurface*)mProperty;
  }

  virtual void SetVisible( bool bVisible = true );
  virtual bool IsVisible();

  int GetVertexIndexAtRAS( double* ras, double* distance );

  int GetVertexIndexAtTarget( double* ras, double* distance, int surface_type = -1 );

  bool GetRASAtVertex       ( int nVertex, double* ras_out, int surface_type = -1 );
  bool GetSurfaceRASAtVertex( int nVertex, double* ras_out );

  int GetVertexAtSurfaceRAS(double* ras, double* distance );

  bool GetTargetAtVertex( int nVertex, double* ras, int surface_type = -1 );

  void GetSurfaceRASAtTarget( double* pos_in, double* ras_out );

  void GetRASAtTarget(double* pos_in, double* ras_out);

  void GetTargetAtSurfaceRAS( double* ras_in, double* pos_out );

  void GetSurfaceRASAtRAS(double* ras_in, double* tkras_out);

  FSSurface* GetSourceSurface()
  {
    return m_surfaceSource;
  }

  void SetPatchFileName( const QString& fn )
  {
    m_sPatchFilename = fn;
  }

  QString GetPatchFileName()
  {
    return m_sPatchFilename;
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

  QList<SurfaceOverlay*> GetOverlays()
  {
    return m_overlays;
  }

  // annotation functions
  int GetNumberOfAnnotations();

  SurfaceAnnotation* GetAnnotation( int n );

  SurfaceAnnotation* GetAnnotation( const QString& name );

  int GetActiveAnnotationIndex();

  SurfaceAnnotation* GetActiveAnnotation();

  void SetActiveAnnotation( int n );

  void SetActiveAnnotation( const QString& name );

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
  void SetActiveLabel(SurfaceLabel* label);

  void DeleteLabel(SurfaceLabel* label);
  void MoveLabelToTop(SurfaceLabel* label);
  void MoveLabelUp(SurfaceLabel* label);
  void MoveLabelDown(SurfaceLabel* label);

  void SetRefVolume(LayerMRI* ref);

  LayerMRI* GetRefVolume()
  {
    return m_volumeRef;
  }

  int GetNumberOfSplines()
  {
    return m_splines.size();
  }

  SurfaceSpline* GetSpline(int n);

  SurfaceSpline* GetActiveSpline()
  {
    return ( m_nActiveSpline >= 0 ? m_splines[m_nActiveSpline] : NULL );
  }

  void DeleteSpline(SurfaceSpline* spline);

  void SetActiveSpline( int n );
  void SetActiveSpline(SurfaceSpline* spline);

  int GetHemisphere();

  void RepositionSurface( LayerMRI* mri, int nVertex, double value, int size, double sigma, int flags = 0 );
  void RepositionSurface( LayerMRI* mri, int nVertex, double* pos, int size, double sigma, int flags = 0 );
  void RepositionSmoothSurface(int nVertex, int size, int n_steps);
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

  int GetCurrentVertex()
  {
    return m_nCurrentVertex;
  }

  int GetLastMark();

  QVector<int> GetAllMarks();

  bool GetCorrelationOverlayDataAtVertex(int nVert, float* output, int nFrames);

  bool IsInflated();

  bool GetActiveLabelCentroidPosition(double* pos);

  bool GetVisibleIn3D()
  {
    return m_bVisibleIn3D;
  }

  void GetSmoothedVertexNormal(int nVertex, double* v_out);

  int GetActiveRGBMap()
  {
    return m_nActiveRGBMap;
  }

  void SetActiveRGBMap(int n);

  QStringList GetRGBMapNames();

  int GetNumberOfRGBMaps()
  {
    return m_rgbMaps.size();
  }

  QString GetMappingSurfaceName()
  {
    return m_sMappingSurfaceName;
  }

  void SetNeighborhoodSize(int nSize);

  QVector<int> GetVertexNeighbors(int vno);

  bool IsContralateralReady()
  {
    return (m_surfaceContralateral != NULL);
  }

  bool IsContralateralPossible();

  LayerSurface* GetContralateralSurface()
  {
    return m_surfaceContralateral;
  }

  int GetContralateralVertex(int vno);

  int GetMouseVertex()
  {
    return m_nMouseVertex;
  }

  void SetMouseVertex(int n)
  {
    m_nMouseVertex = n;
  }

  void AddPathPoint(int vno);

  void RemovePathPoint(int vno);

  void RemoveLastPathPoint();

  void SetActivePath(int n);

  SurfacePath* GetActivePath();

  SurfacePath* GetMadePath(int nPath);

  void DeleteActivePath();

  int FindPathAt(int vno);

  bool IsVertexRipped(int vno);

  bool HasUndoableCut();

  bool FillUncutArea(int vno);

  bool LoadPatch(const QString& filename);

  bool WritePatch(const QString& filename);

  bool FillPath(int nvo, const QVariantMap& options);

  int FillPath(const QVector<int>& verts, const QVariantMap& options);

  void ClearMarks();

  void SaveMarks(const QString& filename);

  SurfaceLabel* CreateNewLabel(const QString& name = "");

  bool LoadParameterization(const QString& filename);

  bool LoadCoordsFromParameterization(const QString &filename);

  void SetSphereFileName(const QString& fn)
  {
    m_sSphereFilename = fn;
  }

  SurfaceAnnotation* CreateNewAnnotation(const QString& ct_file, const QString& name = "");

public slots:
  void SetActiveSurface( int nSurfaceType );
  void UpdateOverlay(bool bAskRedraw = true, bool pre_cached = false);
  void SetLoadAllSurfaces(bool bLoadAll)
  {
    m_bLoadAll = bLoadAll;
  }

  void SetLoadSupSurfaces(const QStringList& names)
  {
    m_listSupFiles = names;
  }

  void SetActiveLabelColor(const QColor& c);
  void SetActiveLabelOutline(bool bOutline);

  void SetActiveAnnotationOutline(bool bOutline);

  void ResetVolumeRef();

  void SetCurrentVertex(int n);

  void UpdateColorMap();

  void RemoveCurrentOverlay();

  void SetVisibleIn3D(bool bVisible);

  void SetHideIn3D(bool bHide)
  {
    SetVisibleIn3D(!bHide);
  }

  void SetMappingSurfaceName(const QString& name)
  {
    if (name.isEmpty())
      m_sMappingSurfaceName = "white";
    else
      m_sMappingSurfaceName = name;
  }

  void AddMappedLabel(LayerROI* label);

  void RemoveMappedLabel(QObject* label_in);

  QVector<int> FindPath(const QVector<int>& seeds);

  void UpdateOverlayLabels()
  {
    UpdateOverlay(true, true);
  }

  void SetContralateralLayer(LayerSurface* layer, LayerSurface* sphere1, LayerSurface* sphere2);
  void ResetContralateralInfo();

  void OnPathCut();

  void OnPathMade();

  void ClearAllCuts();

  void UndoCut();

  QVector<int> FloodFillFromSeed(int seed_vno, const QVariantMap& options = QVariantMap());

  bool IsVertexOnPath(int vno);

  SurfacePath* GetMarks()
  {
    return m_marks;
  }

  vtkActor* GetMainActor();

  bool SaveTransform(const QString& filename);

  void GetCenterOfActor(double* pt);

  bool SavePathAsControlPoints(const QString& fn, bool bMarks = false);

  void SetNoShading(bool b);

Q_SIGNALS:
  void SurfaceAnnotationAdded( SurfaceAnnotation* );
  void SurfaceLabelAdded( SurfaceLabel* );
  void SurfaceLabelDeleted( SurfaceLabel* );
  void SurfaceSplineAdded( SurfaceSpline* );
  void SurfaceSplineDeleted( SurfaceSpline* );
  void SurfaceRGBAdded();
  void SurfaceOverlayAdded( SurfaceOverlay* );
  void SurfaceOverlyDataUpdated();
  void SurfaceCurvatureLoaded();
  void SurfaceVectorLoaded();
  void ActiveSurfaceChanged( int n );
  void ActiveOverlayChanged( int n );
  void ActiveAnnotationChanged( int n );
  void ActiveLabelChanged( int n );
  void ActiveSplineChanged( int n );
  void CurrentVertexChanged(int n);
  void RGBMapChanged();
  void FlattenedPatchLoaded();

protected slots:
  void UpdateOpacity();
  void UpdateEdgeThickness();
  void UpdateVectorPointSize();
  void UpdateRenderMode();
  void UpdateVertexRender();
  void UpdateMeshRender();
  void UpdateActorPositions();
  void UpdateROIPosition(double dx, double dy, double dz);
  void UpdateVectorActor2D();
  void OnSlicePositionChanged3D();
  void SetHighlightedLabelOnAnnotation(int n);

protected:
  void InitializeData();
  void InitializeSurface();
  void InitializeActors();
  void InitializeLabel(SurfaceLabel* label);
  void MapLabels( unsigned char* data, int nVertexCount );
  void EditPathPoint(int vno, bool remove = false);
  void PushMarksToPath();
  void DoSlicePositionChanged(int nPlane, bool bUpdatePosOnly = false);

  virtual void OnSlicePositionChanged( int nPlane );
  virtual void OnSetDisplayInNeurologicalView();

  // Pipeline ------------------------------------------------------------
  vtkSmartPointer<vtkPlane>     mReslicePlane[3];
  vtkSmartPointer<vtkImageMapToColors>  mColorMap[3];

  FSSurface*   m_surfaceSource;
  bool    m_bResampleToRAS;
  LayerMRI*   m_volumeRef;

  QString   m_sPatchFilename;
  QString   m_sVectorFilename;
  QString   m_sTargetFilename;
  QString   m_sSphereFilename;

  vtkSmartPointer<vtkActor>   m_sliceActor2D[3];
  vtkSmartPointer<vtkActor>   m_sliceActor3D[3];
  vtkSmartPointer<vtkActor>   m_vectorActor2D[3];

  vtkSmartPointer<vtkActor>   m_mainActor;
  vtkSmartPointer<vtkActor>   m_vectorActor;
  vtkSmartPointer<vtkActor>   m_vertexActor;
  vtkSmartPointer<vtkActor>   m_vertexActor2D[3];
  vtkSmartPointer<vtkActor>   m_wireframeActor;

  vtkSmartPointer<vtkCutter>  m_cutter[3];
  vtkSmartPointer<vtkBox>     m_box[3];
  vtkSmartPointer<vtkPolyData>     m_vertexPoly2D[3];

  QList<SurfaceOverlay*>    m_overlays;
  int         m_nActiveOverlay;

  QList<SurfaceAnnotation*> m_annotations;
  int         m_nActiveAnnotation;

  QList<SurfaceLabel*>      m_labels;
  int         m_nActiveLabel;

  QList<SurfacePath*> m_paths;
  int         m_nActivePath;
  SurfacePath*  m_marks;

  QList<RGBMap>             m_rgbMaps;
  int         m_nActiveRGBMap;

  int         m_nCurrentVertex;
  int         m_nMouseVertex;

  SurfaceROI*           m_roi;
  QList<SurfaceSpline*> m_splines;
  int         m_nActiveSpline;

  bool        m_bUndoable;
  bool        m_bVector2DPendingUpdate;
  bool        m_bLoadAll;

  bool        m_bVisibleIn3D;

  QString     m_sMappingSurfaceName;

  QList<LayerROI*>  m_mappedLabels;

  QStringList m_listSupFiles;

  unsigned char*    m_nColorDataCache;

  LayerSurface*     m_surfaceContralateral;
  LayerSurface*     m_surfaceSphere1;
  LayerSurface*     m_surfaceSphere2;
};

#endif


