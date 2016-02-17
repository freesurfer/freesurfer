/**
 * @file  LayerSurface.cpp
 * @brief Layer data object for MRI surface.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: zkaufman $
 *    $Date: 2016/02/17 20:36:46 $
 *    $Revision: 1.113 $
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

#include <QSharedPointer>
#include "MainWindow.h"
#include "LayerSurface.h"
#include "vtkRenderer.h"
#include "vtkImageActor.h"
#include "vtkPolyDataMapper.h"
#include "vtkSmartPointer.h"
#include "vtkMatrix4x4.h"
#include "vtkImageMapToColors.h"
#include "vtkImageActor.h"
#include "vtkActor.h"
#include "vtkLODActor.h"
#include "vtkRGBAColorTransferFunction.h"
#include "vtkLookupTable.h"
#include "vtkProperty.h"
#include "vtkImageReslice.h"
#include "vtkFreesurferLookupTable.h"
#include "vtkPlane.h"
#include "vtkCutter.h"
#include "vtkDecimatePro.h"
#include "vtkMapperCollection.h"
#include "vtkTubeFilter.h"
#include "vtkPointData.h"
#include "vtkTubeFilter.h"
#include "FSSurface.h"
#include "LayerMRI.h"
#include "SurfaceAnnotation.h"
#include "SurfaceLabel.h"
#include "LayerPropertySurface.h"
#include "SurfaceROI.h"
#include <QFileInfo>
#include <QDir>
#include <QDebug>
#include "SurfaceOverlayProperty.h"
#include "MyUtils.h"
#include "SurfaceOverlay.h"
#include "SurfaceSpline.h"
#include "vtkMaskPoints.h"
#include "vtkExtractPolyDataGeometry.h"
#include "vtkBox.h"

LayerSurface::LayerSurface( LayerMRI* ref, QObject* parent ) : LayerEditable( parent ),
  m_surfaceSource( NULL ),
  m_bResampleToRAS( true ),
  m_volumeRef( ref ),
  m_nActiveOverlay( -1 ),
  m_nActiveAnnotation( -1 ),
  m_nActiveLabel( -1 ),
  m_bUndoable( false ),
  m_bVector2DPendingUpdate( true ),
  m_bLoadAll(false),
  m_spline(NULL),
  m_nCurrentVertex(-1)
{
  m_strTypeNames.push_back( "Surface" );
  m_sPrimaryType = "Surface";

  // create property before actors!
  mProperty = new LayerPropertySurface( this );

  for ( int i = 0; i < 3; i++ )
  {
    // m_nSliceNumber[i] = 0;
    m_sliceActor2D[i] = vtkSmartPointer<vtkActor>::New();
    m_sliceActor3D[i] = vtkSmartPointer<vtkActor>::New();
    m_vectorActor2D[i] = vtkSmartPointer<vtkActor>::New();
    m_vertexActor2D[i] = vtkSmartPointer<vtkActor>::New();
    m_vertexActor2D[i]->GetProperty()->SetRepresentationToPoints();
    m_vertexActor2D[i]->VisibilityOff();
  }

  m_mainActor = vtkSmartPointer<vtkActor>::New();
  m_mainActor->GetProperty()->SetEdgeColor( 0.75, 0.75, 0.75 );

  m_vectorActor = vtkSmartPointer<vtkActor>::New();
  m_vectorActor->GetProperty()->SetColor( GetProperty()->GetVectorColor() );
  m_vectorActor->GetProperty()->SetPointSize( GetProperty()->GetVectorPointSize() );
  m_vectorActor->PickableOff();

  for ( int i = 0; i < 3; i++ )
  {
    m_vectorActor2D[i]->GetProperty()->SetColor( GetProperty()->GetVectorColor() );
    m_vectorActor2D[i]->GetProperty()->SetPointSize( GetProperty()->GetVectorPointSize() );
    m_vectorActor2D[i]->PickableOff();
  }

  m_vertexActor = vtkSmartPointer<vtkActor>::New();
  m_vertexActor->GetProperty()->SetRepresentationToPoints();
  m_vertexActor->VisibilityOff();

  m_wireframeActor = vtkSmartPointer<vtkActor>::New();
  m_wireframeActor->VisibilityOff();

  m_roi = new SurfaceROI(this);

  m_spline = new SurfaceSpline(this);
  connect(m_spline, SIGNAL(SplineChanged()), this, SIGNAL(ActorChanged()));

  LayerPropertySurface* p = GetProperty();
  connect( p, SIGNAL(ColorMapChanged()), this, SLOT(UpdateColorMap()) );
  connect( p, SIGNAL(OpacityChanged(double)), this, SLOT(UpdateOpacity()) );
  connect( p, SIGNAL(EdgeThicknessChanged(int)), this, SLOT(UpdateEdgeThickness()) );
  connect( p, SIGNAL(VectorPointSizeChanged(int)), this, SLOT(UpdateVectorPointSize()) );
  connect( p, SIGNAL(RenderModeChanged(int)), this, SLOT(UpdateRenderMode()) );
  connect( p, SIGNAL(VertexRenderChanged()), this, SLOT(UpdateVertexRender()) );
  connect( p, SIGNAL(MeshRenderChanged()), this, SLOT(UpdateMeshRender()) );
  connect( p, SIGNAL(PositionChanged()), this, SLOT(UpdateActorPositions()) );
  connect( p, SIGNAL(PositionChanged(double, double, double)),
           this, SLOT(UpdateROIPosition(double, double, double)));
  connect( p, SIGNAL(OverlayChanged()), this, SLOT(UpdateOverlay()));
  connect(this, SIGNAL(ActiveLabelChanged(int)), p, SIGNAL(PropertyChanged()));

  if (m_volumeRef)
    connect( m_volumeRef, SIGNAL(destroyed()), this, SLOT(ResetVolumeRef()), Qt::UniqueConnection);
}

LayerSurface::~LayerSurface()
{
  if ( m_surfaceSource )
  {
    delete m_surfaceSource;
  }

  for ( int i = 0; i < m_overlays.size(); i++ )
  {
    delete m_overlays[i];
  }

  for ( int i = 0; i < m_annotations.size(); i++ )
  {
    delete m_annotations[i];
  }

  for ( int i = 0; i < m_labels.size(); i++ )
  {
    delete m_labels[i];
  }
}

void LayerSurface::SetRefVolume(LayerMRI *ref)
{
  m_volumeRef = ref;
  if (m_volumeRef)
    connect( m_volumeRef, SIGNAL(destroyed()), this, SLOT(ResetVolumeRef()), Qt::UniqueConnection);
}

bool LayerSurface::LoadSurfaceFromFile()
{
  if ( m_surfaceSource )
  {
    delete m_surfaceSource;
  }

  m_surfaceSource = new FSSurface( m_volumeRef ? m_volumeRef->GetSourceVolume() : NULL );
  if ( !m_surfaceSource->MRISRead( m_sFilename,
                                   m_sVectorFilename,
                                   m_sPatchFilename,
                                   m_sTargetFilename,
                                   m_listSupFiles )
     )
  {
    return false;
  }

  InitializeData();

  return true;
}

void LayerSurface::InitializeData()
{
  ParseSubjectName(m_sFilename);
  InitializeSurface();
  InitializeActors();

  GetProperty()->SetSurfaceSource( m_surfaceSource );

  if (IsInflated())
  {
    double pos[3] = { -45, 0, 0 };
    if (GetHemisphere() == 1)
      pos[0] = -pos[0];
    GetProperty()->SetPosition(pos);
    GetProperty()->SetEdgeThickness(0);
  }
}

bool LayerSurface::CreateFromMRIS(void *mris_ptr)
{
  MRIS* mris = (MRIS*)mris_ptr;
  m_surfaceSource = new FSSurface( m_volumeRef ? m_volumeRef->GetSourceVolume() : NULL );
  if ( !m_surfaceSource->CreateFromMRIS(mris) )
  {
    return false;
  }
  SetFileName(mris->fname);
  InitializeData();
  return true;
}

bool LayerSurface::SaveSurface( const QString& filename )
{
  if ( !m_surfaceSource->MRISWrite( filename ) )
  {
    cerr << "MRISWrite failed: Unable to write to " << qPrintable(filename) << ".\n";
    return false;
  }
  else
  {
    ResetModified();
    return true;
  }
}

bool LayerSurface::SaveSurface( )
{
  if ( m_sFilename.size() == 0 )
  {
    cerr << "No filename provided to save surface.\n";
    return false;
  }

  return SaveSurface( m_sFilename.toAscii().data() );
}

bool LayerSurface::LoadVectorFromFile( )
{
  if ( m_sVectorFilename.size() == 0 || !m_surfaceSource->MRISReadVectors( m_sVectorFilename ) )
  {
    return false;
  }

  UpdateVectorActor2D();
  emit Modified();
  emit SurfaceVectorLoaded();
  emit ActorUpdated();

  return true;
}

void LayerSurface::UpdateVectorActor2D()
{
  if ( m_surfaceSource && m_surfaceSource->GetActiveVector() >= 0 )
  {
    for ( int i = 0; i < 3; i++ )
    {
      vtkPolyDataMapper* mapper = vtkPolyDataMapper::SafeDownCast( m_sliceActor2D[i]->GetMapper() );
      if ( mapper )
      {
        mapper->Update();
      }

      m_surfaceSource->UpdateVector2D( i,
                                       m_dSlicePosition[i],
                                       ( mapper ? mapper->GetInput() : NULL )
                                     );
    }
  }
}

bool LayerSurface::LoadCurvatureFromFile( const QString& filename )
{
  QString fn = filename;
  fn.replace("~", QDir::homePath());
  if (!QFile::exists(fn))
  {
    fn = QFileInfo(QFileInfo(m_sFilename).dir(), filename).absoluteFilePath();
  }

  if ( !m_surfaceSource->LoadCurvature( fn ) )
  {
    return false;
  }

  UpdateOverlay(false);
  emit Modified();
  emit SurfaceCurvatureLoaded();
  emit ActorUpdated();

  return true;
}

bool LayerSurface::LoadOverlayFromFile(const QString &filename, const QString& fn_reg, bool bCorrelation, bool bSeconfHalfData)
{
  QString fn = filename;
  fn.replace("~", QDir::homePath());
  fn = QFileInfo(fn).absoluteFilePath();
  if (!QFile::exists(fn))
  {
    fn = QFileInfo(QFileInfo(m_sFilename).dir(), filename).absoluteFilePath();
  }

  if (bCorrelation)
  {
    return LoadCorrelationFromFile(fn);
  }
  else
  {
    QString fullpath = fn_reg;
    if (!fn_reg.isEmpty())
    {
      fullpath.replace("~", QDir::homePath());
      fullpath = QFileInfo(fullpath).absoluteFilePath();
    }
    return LoadGenericOverlayFromFile(fn, fullpath, bSeconfHalfData);
  }
}

bool LayerSurface::LoadGenericOverlayFromFile( const QString& filename, const QString& fn_reg, bool bSecondHalfData )
{
  float* data = NULL;
  int nframes, nvertices;
  if ( !m_surfaceSource->LoadOverlay( filename, fn_reg, &data, &nvertices, &nframes, bSecondHalfData ) )
  {
    return false;
  }

  // create overlay
  SurfaceOverlay* overlay = new SurfaceOverlay( this );
  overlay->InitializeData(data, nvertices, nframes);
  overlay->SetName( QFileInfo(filename).fileName() );
  overlay->SetFileName( filename );

  m_overlays.push_back( overlay );
  SetActiveOverlay( m_overlays.size() - 1 );

  emit Modified();
  emit SurfaceOverlayAdded( overlay );
  connect(overlay, SIGNAL(DataUpdated()), this, SIGNAL(SurfaceOverlyDataUpdated()), Qt::UniqueConnection);
  return true;
}

bool LayerSurface::LoadCorrelationFromFile( const QString& filename )
{
  // create overlay
  SurfaceOverlay* overlay = new SurfaceOverlay( this );
  overlay->SetName( QFileInfo(filename).fileName() );
  overlay->SetFileName( filename );
  if ( !overlay->LoadCorrelationData( filename ) )
  {
    return false;
  }

  m_overlays.push_back( overlay );
  SetActiveOverlay( m_overlays.size() - 1 );

  emit Modified();
  emit SurfaceOverlayAdded( overlay );
  connect(overlay, SIGNAL(DataUpdated()), this, SIGNAL(SurfaceOverlyDataUpdated()), Qt::UniqueConnection);
  return true;
}

void LayerSurface::CopyCorrelationOverlay(LayerSurface *surf)
{
  SurfaceOverlay* src = surf->GetCorrelationOverlay();
  if (src)
  {
    SurfaceOverlay* overlay = new SurfaceOverlay( this );
    overlay->SetName(src->GetName());
    overlay->CopyCorrelationData(src);
    m_overlays.push_back( overlay );
    SetActiveOverlay( m_overlays.size() - 1 );

    emit Modified();
    emit SurfaceOverlayAdded( overlay );
    connect(overlay, SIGNAL(DataUpdated()), this, SIGNAL(SurfaceOverlyDataUpdated()), Qt::UniqueConnection);
  }
}

bool LayerSurface::LoadAnnotationFromFile( const QString& filename )
{
  QString fn = filename;
  if (fn.left(2) == "~/")
    fn.replace("~", QDir::homePath());
  QFileInfo fi(fn);

  // create annotation
  SurfaceAnnotation* annot = new SurfaceAnnotation( this );
  bool ret = annot->LoadAnnotation( fn );
  if ( !ret )
  {
    delete annot;
    return false;
  }

  QString name;
  if ( fi.suffix() == ".annot" )
  {
    name = fi.completeBaseName();
  }
  else
  {
    name = fi.fileName();
  }

  QStringList names;
  for (int i = 0; i < m_annotations.size(); i++)
    names << m_annotations[i]->GetName();

  QString basename = name;
  int n = 0;
  while (names.contains(name))
  {
    n++;
    name = QString("%1_%2").arg(basename).arg(n);
  }

  annot->SetName(name);
  m_annotations.push_back( annot );

  SetActiveAnnotation( m_annotations.size() - 1 );

  emit Modified();
  emit SurfaceAnnotationAdded( annot );
  emit ActorUpdated();
  return true;
}

bool LayerSurface::LoadLabelFromFile( const QString& filename )
{
  // create label
  SurfaceLabel* label = new SurfaceLabel( this );

  QString fn = filename;
  if (fn.left(2) == "~/")
    fn.replace("~", QDir::homePath());
  QFileInfo fi(fn);
  bool ret = label->LoadLabel(fn);
  if ( !ret )
  {
    delete label;
    return false;
  }

  if ( fi.suffix() == ".label" )
  {
    label->SetName( fi.completeBaseName() );
  }
  else
  {
    label->SetName( fi.fileName() );
  }

  m_labels.insert(0, label);
  connect(label, SIGNAL(SurfaceLabelChanged()), this, SLOT(UpdateColorMap()));
  connect(label, SIGNAL(SurfaceLabelChanged()), GetProperty(), SIGNAL(PropertyChanged()));
  connect(label, SIGNAL(SurfaceLabelVisibilityChanged()), this, SLOT(UpdateColorMap()));

  SetActiveLabel( 0 );

  UpdateOverlay(false);

  emit Modified();
  emit SurfaceLabelAdded( label );
  emit ActorChanged();
  return true;
}

bool LayerSurface::LoadSplineFromFile(const QString &filename)
{
  if (!m_spline->Load(filename))
  {
    return false;
  }

  return true;
}

void LayerSurface::InitializeSurface()
{
  if ( m_surfaceSource == NULL )
  {
    return;
  }

  FSSurface* source = m_surfaceSource;

  float RASBounds[6];
  source->GetBounds( RASBounds );
  m_dWorldOrigin[0] = RASBounds[0];
  m_dWorldOrigin[1] = RASBounds[2];
  m_dWorldOrigin[2] = RASBounds[4];

  m_dWorldSize[0] = RASBounds[1] - RASBounds[0];
  m_dWorldSize[1] = RASBounds[3] - RASBounds[2];
  m_dWorldSize[2] = RASBounds[5] - RASBounds[4];
}

void LayerSurface::InitializeActors()
{
  if ( m_surfaceSource == NULL )
  {
    return;
  }

  // main surface actor
  vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  mapper->SetInput( m_surfaceSource->GetPolyData() );
  m_mainActor->SetMapper( mapper );
  mapper->Update();

  // vector actor
  mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  vtkSmartPointer<vtkTubeFilter> tube = vtkSmartPointer<vtkTubeFilter>::New();
  tube->SetInput(m_surfaceSource->GetVectorPolyData());
  tube->SetNumberOfSides(8);
  tube->SetRadius(0.04);
  tube->CappingOn();
  mapper->SetInput( tube->GetOutput() );
  m_vectorActor->SetMapper( mapper );
//  mapper->Update();

  for ( int i = 0; i < 3; i++ )
  {
    mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInput(  m_surfaceSource->GetVector2DPolyData( i ) );
    m_vectorActor2D[i]->SetMapper( mapper );
  }

  // vertex actor
  mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  mapper->SetInput(  m_surfaceSource->GetVertexPolyData() );
  m_vertexActor->SetMapper( mapper );

  // wireframe actor
  mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  mapper->SetInput(  m_surfaceSource->GetWireframePolyData() );
  m_wireframeActor->SetMapper( mapper );
  mapper->Update();

  for ( int i = 0; i < 3; i++ )
  {
    // The reslice object just takes a slice out of the volume.
    //
    mReslicePlane[i] = vtkSmartPointer<vtkPlane>::New();
    mReslicePlane[i]->SetOrigin( 0, 0, 0 );
    mReslicePlane[i]->SetNormal( (i==0), (i==1), (i==2) );

    double bounds[6];
    m_surfaceSource->GetPolyData()->GetBounds(bounds);
    {
      bounds[i*2] = (bounds[i*2]+bounds[i*2+1])/2-0;
      bounds[i*2+1] = bounds[i*2]+20;
    }
    m_box[i] = vtkSmartPointer<vtkBox>::New();
    m_box[i]->SetBounds(bounds);
    vtkSmartPointer<vtkExtractPolyDataGeometry> extract = vtkSmartPointer<vtkExtractPolyDataGeometry>::New();
    extract->SetInput(m_surfaceSource->GetPolyData());
    extract->ExtractInsideOn();
    extract->ExtractBoundaryCellsOn();
    extract->SetImplicitFunction(m_box[i]);

    m_cutter[i] =
      vtkSmartPointer<vtkCutter>::New();
    m_cutter[i]->SetInputConnection( extract->GetOutputPort() );
    m_cutter[i]->SetCutFunction( mReslicePlane[i] );
    m_cutter[i]->GenerateCutScalarsOff();

    //
    // Mappers for the lines.
    //
    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputConnection( m_cutter[i]->GetOutputPort() );
//   mapper->SetInputConnection( 1, cutter->GetOutputPort( 1 ) );
    vtkSmartPointer<vtkPolyDataMapper> mapper2 = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper2->SetInputConnection( m_cutter[i]->GetOutputPort() );
    //
    // Actors in the scene, drawing the mapped lines.
    //
    m_sliceActor2D[i]->SetMapper( mapper );
//  m_sliceActor2D[i]->SetBackfaceProperty( m_sliceActor2D[i]->MakeProperty() );
//  m_sliceActor2D[i]->GetBackfaceProperty()->BackfaceCullingOff();
    m_sliceActor2D[i]->SetProperty( m_sliceActor2D[i]->MakeProperty() );
    m_sliceActor2D[i]->GetProperty()->SetInterpolationToFlat();
    m_sliceActor2D[i]->GetProperty()->SetLineWidth( GetProperty()->GetEdgeThickness() );

    m_sliceActor3D[i]->SetMapper( mapper2 );
//  m_sliceActor3D[i]->SetBackfaceProperty( m_sliceActor3D[i]->MakeProperty() );
//  m_sliceActor3D[i]->GetBackfaceProperty()->BackfaceCullingOff();
    m_sliceActor3D[i]->SetProperty( m_sliceActor3D[i]->MakeProperty() );
    m_sliceActor3D[i]->GetProperty()->SetLineWidth( GetProperty()->GetEdgeThickness() );
//    m_sliceActor3D[i]->GetProperty()->SetInterpolationToFlat();

    vtkSmartPointer<vtkPolyDataMapper> mapper3 = vtkSmartPointer<vtkPolyDataMapper>::New();
//    vtkSmartPointer<vtkMaskPoints> pts = vtkSmartPointer<vtkMaskPoints>::New();
//    pts->GenerateVerticesOn();
//    pts->SetOnRatio(1);
//    pts->SetInputConnection(m_cutter[i]->GetOutputPort());
    mapper3->SetInputConnection( m_cutter[i]->GetOutputPort() );
    mapper3->ScalarVisibilityOff();
    m_vertexActor2D[i]->SetMapper(mapper3);
    m_vertexActor2D[i]->SetProperty( m_vertexActor2D[i]->MakeProperty() );
    m_vertexActor2D[i]->GetProperty()->SetPointSize(3);
    m_vertexActor2D[i]->GetProperty()->SetInterpolationToFlat();

    // Set ourselves up.
    this->OnSlicePositionChanged( i );
  }

  this->UpdateOpacity();
  this->UpdateColorMap();
  this->UpdateVectorActor2D();
}

void LayerSurface::UpdateOpacity()
{
  for ( int i = 0; i < 3; i++ )
  {
    // m_sliceActor2D[i]->GetProperty()->SetOpacity( mProperty->GetOpacity() );
    // m_sliceActor3D[i]->SetOpacity( mProperty->GetOpacity() );
  }
  m_mainActor->GetProperty()->SetOpacity( GetProperty()->GetOpacity() );
  emit ActorUpdated();
}

void LayerSurface::UpdateEdgeThickness()
{
  for ( int i = 0; i < 3; i++ )
  {
    m_sliceActor2D[i]->GetProperty()->SetLineWidth( GetProperty()->GetEdgeThickness() );
    m_sliceActor3D[i]->GetProperty()->SetLineWidth( GetProperty()->GetEdgeThickness() );
  }

  if ( GetProperty()->GetEdgeThickness() == 0 )
  {
    if ( IsVisible() )
    {
      for ( int i = 0; i < 3; i++ )
      {
        m_sliceActor2D[i]->SetVisibility( 0 );
        m_sliceActor3D[i]->SetVisibility( 0 );
      }
    }
  }
  else if ( IsVisible() && !m_sliceActor2D[0]->GetVisibility() )
  {
    for ( int i = 0; i < 3; i++ )
    {
      m_sliceActor2D[i]->SetVisibility( 1 );
      m_sliceActor3D[i]->SetVisibility( 1 );
    }
  }
  emit ActorUpdated();
}

void LayerSurface::UpdateVectorPointSize()
{
  for ( int i = 0; i < 3; i++ )
  {
    m_vectorActor2D[i]->GetProperty()->SetPointSize( GetProperty()->GetVectorPointSize() );
  }
  m_vectorActor->GetProperty()->SetPointSize( GetProperty()->GetVectorPointSize() );
  emit ActorUpdated();
}

void LayerSurface::UpdateColorMap()
{
  if ( m_surfaceSource == NULL )
  {
    return;
  }

  for ( int i = 0; i < 3; i++ )
  {
    m_sliceActor2D[i]->GetProperty()->SetColor( GetProperty()->GetEdgeColor() );
    m_sliceActor3D[i]->GetProperty()->SetColor( GetProperty()->GetEdgeColor() );
    m_sliceActor2D[i]->GetMapper()->ScalarVisibilityOff();
    m_sliceActor3D[i]->GetMapper()->ScalarVisibilityOff();
    m_vectorActor2D[i]->GetProperty()->SetColor( GetProperty()->GetVectorColor() );
  }

  m_mainActor->GetProperty()->SetColor( GetProperty()->GetBinaryColor() );
  m_wireframeActor->GetProperty()->SetColor( GetProperty()->GetBinaryColor() );
  m_vectorActor->GetProperty()->SetColor( GetProperty()->GetVectorColor() );
  if ( m_surfaceSource->IsCurvatureLoaded() )
  {
    if ( GetProperty()->GetCurvatureLUT() != m_mainActor->GetMapper()->GetLookupTable() )
    {
      m_mainActor->GetMapper()->SetLookupTable( GetProperty()->GetCurvatureLUT() );

    }

    if ( GetProperty()->GetMeshColorMap() == LayerPropertySurface::MC_Surface )
    {
      m_wireframeActor->GetMapper()->SetLookupTable( GetProperty()->GetCurvatureLUT() );
    }
    else if ( GetProperty()->GetMeshColorMap() == LayerPropertySurface::MC_Curvature )
    {
      UpdateMeshRender();
    }
    /*  vtkSmartPointer<vtkMapperCollection> mc = m_mainActor->GetLODMappers();
      mc->InitTraversal();
      vtkMapper* mapper = NULL;
      while ( ( mapper = mc->GetNextItem() ) != NULL )
      {
    mapper->SetLookupTable( GetProperty()->GetCurvatureLUT() );
      } */
  }

  UpdateOverlay(false);
  emit ActorUpdated();
}

void LayerSurface::Append2DProps( vtkRenderer* renderer, int nPlane )
{
  renderer->AddViewProp( m_sliceActor2D[nPlane] );
  renderer->AddViewProp( m_vectorActor2D[nPlane] );
  renderer->AddViewProp( m_vertexActor2D[nPlane]);
  if (m_spline)
    m_spline->AppendProp2D(renderer, nPlane);
}

void LayerSurface::Append3DProps( vtkRenderer* renderer, bool* bSliceVisibility )
{
  for ( int i = 0; i < 3; i++ )
  {
    if (bSliceVisibility[i])
    {
      renderer->AddViewProp( m_sliceActor3D[i] );
    }
  }

  renderer->AddViewProp( m_mainActor );
  renderer->AddViewProp( m_vectorActor );
  renderer->AddViewProp( m_vertexActor );
  renderer->AddViewProp( m_wireframeActor );

  if (m_spline)
    m_spline->AppendProp3D(renderer);

  m_roi->AppendProps(renderer);
}

/*
void LayerSurface::SetSliceNumber( int* sliceNumber )
{
 if ( sliceNumber[0] != m_nSliceNumber[0] || sliceNumber[1] != m_nSliceNumber[1] ||
   sliceNumber[2] != m_nSliceNumber[2] )
 {
  m_nSliceNumber[0] = sliceNumber[0];
  m_nSliceNumber[1] = sliceNumber[1];
  m_nSliceNumber[2] = sliceNumber[2];


 }
}
*/

void LayerSurface::SetSlicePositionToWorldCenter()
{
  if ( m_surfaceSource == NULL )
  {
    return;
  }

  // Get some values from the MRI.
  double pos[3];
  for ( int i = 0; i < 3; i++ )
  {
    pos[i] = ((int)( m_dWorldSize[i]/2/m_dWorldVoxelSize[i] ) + 0.0 ) * m_dWorldVoxelSize[i] + m_dWorldOrigin[i];
  }

  SetSlicePosition( pos );
}

void LayerSurface::OnSlicePositionChanged( int nPlane )
{
  if ( m_surfaceSource == NULL )
  {
    return;
  }

  double* pos = GetProperty()->GetPosition();
  double bounds[6];
  m_surfaceSource->GetPolyData()->GetBounds(bounds);
  switch ( nPlane )
  {
  case 0:
    mReslicePlane[0]->SetOrigin( m_dSlicePosition[0]-pos[0], 0, 0  );
    m_sliceActor2D[0]->SetPosition( 0.1, pos[1], pos[2] );
    m_vertexActor2D[0]->SetPosition( 0.1, pos[1], pos[2] );
    m_vectorActor2D[0]->SetPosition( 1.0, pos[1], pos[2] );
    break;
  case 1:
    mReslicePlane[1]->SetOrigin( 0, m_dSlicePosition[1]-pos[1], 0 );
    m_sliceActor2D[1]->SetPosition( pos[0], 0.1, pos[2] );
    m_vertexActor2D[1]->SetPosition( pos[0], 0.1, pos[2] );
    m_vectorActor2D[1]->SetPosition( pos[0], 1.0, pos[2] );
    break;
  case 2:
    mReslicePlane[2]->SetOrigin( 0, 0, m_dSlicePosition[2]-pos[2]  );
    m_sliceActor2D[2]->SetPosition( pos[0], pos[1], -0.1 );
    m_vertexActor2D[2]->SetPosition( pos[0], pos[1], -0.1 );
    m_vectorActor2D[2]->SetPosition( pos[0], pos[1], -1.0 );
    break;
  }
  double dLen = m_surfaceSource->GetMaxSegmentLength();
  bounds[nPlane*2] = m_dSlicePosition[nPlane]-dLen/2;
  bounds[nPlane*2+1] = bounds[nPlane*2]+dLen/2;
  m_box[nPlane]->SetBounds(bounds);

  // update mapper so the polydata is current
  if ( IsVisible() && GetActiveVector() >= 0 )
  {
    vtkPolyDataMapper* mapper = vtkPolyDataMapper::SafeDownCast( m_sliceActor2D[nPlane]->GetMapper() );
    if ( mapper )
    {
      mapper->Update();
    }
    m_surfaceSource->UpdateVector2D( nPlane, m_dSlicePosition[nPlane], ( mapper ? mapper->GetInput() : NULL ) );
  }
  else
  {
    m_bVector2DPendingUpdate = true;
  }

  if (m_spline && MainWindow::GetMainWindow()->GetSplinePicking())
  {
    int nVertex = this->GetVertexIndexAtTarget(m_dSlicePosition, NULL);
    m_spline->SetActiveVertex(nVertex);
  }

  for (int i = 0; i < m_overlays.size(); i++)
  {
    m_overlays[i]->UpdateCorrelationCoefficient(m_dSlicePosition);
  }
}

void LayerSurface::SetVisible( bool bVisible )
{
  int nSliceVisibility = ( ( bVisible && GetProperty()->GetEdgeThickness() > 0 ) ? 1 : 0 );
  for ( int i = 0; i < 3; i++ )
  {
    m_sliceActor2D[i] ->SetVisibility( nSliceVisibility );
    m_sliceActor3D[i] ->SetVisibility( nSliceVisibility );
    m_vertexActor2D[i]->SetVisibility( bVisible && GetProperty()->GetShowVertices() );
  }

  m_mainActor->SetVisibility( bVisible && GetProperty()->GetSurfaceRenderMode() != LayerPropertySurface::SM_Wireframe );
  m_wireframeActor->SetVisibility( bVisible && GetProperty()->GetSurfaceRenderMode() != LayerPropertySurface::SM_Surface );
  m_vertexActor->SetVisibility( bVisible && GetProperty()->GetShowVertices() );

  int nVectorVisibility = ( ( bVisible && m_surfaceSource && m_surfaceSource->GetActiveVector() >= 0 )? 1 : 0 );
  m_vectorActor->SetVisibility( nVectorVisibility );
  if ( nVectorVisibility && m_bVector2DPendingUpdate )
  {
    UpdateVectorActor2D();
    m_bVector2DPendingUpdate = false;
  }
  for ( int i = 0; i < 3; i++ )
  {
    m_vectorActor2D[i]->SetVisibility( nVectorVisibility );
  }

  m_spline->SetVisible(bVisible);

  LayerEditable::SetVisible(bVisible);
}

bool LayerSurface::IsVisible()
{
  if ( GetProperty()->GetSurfaceRenderMode() == LayerPropertySurface::SM_Wireframe )
  {
    return m_wireframeActor->GetVisibility() > 0;
  }
  else
  {
    return m_mainActor->GetVisibility() > 0;
  }
}

bool LayerSurface::HasProp( vtkProp* prop )
{
  for ( int i = 0; i < 3; i++ )
  {
    if ( m_sliceActor2D[i].GetPointer() == prop || m_sliceActor3D[i].GetPointer() == prop )
    {
      return true;
    }
  }
  return (m_mainActor.GetPointer() == prop || m_wireframeActor.GetPointer() == prop);
}

int LayerSurface::GetVertexIndexAtRAS( double* ras, double* distance )
{
  if ( m_surfaceSource == NULL )
  {
    return -1;
  }

  return m_surfaceSource->FindVertexAtRAS( ras, distance );
}

int LayerSurface::GetVertexIndexAtTarget( double* pos, double* distance )
{
  if ( m_surfaceSource == NULL )
  {
    return -1;
  }

  double pos_o[3];
  double* offset = GetProperty()->GetPosition();
  for ( int i = 0; i < 3; i++ )
  {
    pos_o[i] = pos[i] - offset[i];
  }
  /*
  if ( m_volumeRef )
  {
    double realRas[3];
    m_volumeRef->TargetToRAS( pos_o, realRas );
    return m_surfaceSource->FindVertexAtRAS( realRas, distance );
  }
  else
  {
    return m_surfaceSource->FindVertexAtRAS( pos_o, distance );
  }
  */
  double realRas[3];
  m_surfaceSource->ConvertTargetToRAS( pos_o, realRas );
  return m_surfaceSource->FindVertexAtRAS( realRas, distance );
}

bool LayerSurface::GetRASAtVertex( int nVertex, double* ras )
{
  if ( m_surfaceSource == NULL )
  {
    return false;
  }

  return m_surfaceSource->GetRASAtVertex( nVertex, ras );
}

void LayerSurface::GetSurfaceRASAtTarget( double* pos_in, double* ras_out )
{
  if ( m_surfaceSource == NULL )
  {
    return;
  }

  double pos_o[3];
  double* offset = GetProperty()->GetPosition();
  for ( int i = 0; i < 3; i++ )
  {
    pos_o[i] = pos_in[i] - offset[i];
  }
  /*
  if ( m_volumeRef )
  {
    m_volumeRef->TargetToRAS( pos_o, pos_o );
  }
  */
  m_surfaceSource->ConvertTargetToRAS( pos_o, pos_o );
  m_surfaceSource->ConvertRASToSurface( pos_o, ras_out );
}

void LayerSurface::GetTargetAtSurfaceRAS( double* ras_in, double* pos_out )
{
  if ( m_surfaceSource == NULL )
  {
    return;
  }

  m_surfaceSource->ConvertSurfaceToRAS( ras_in, pos_out );
  /*
  if ( m_volumeRef )
  {
    m_volumeRef->RASToTarget( pos_out, pos_out );
  }
  */
  m_surfaceSource->ConvertRASToTarget(pos_out, pos_out);
}

bool LayerSurface::GetSurfaceRASAtVertex( int nVertex, double* ras )
{
  if ( m_surfaceSource == NULL )
  {
    return false;
  }

  return m_surfaceSource->GetSurfaceRASAtVertex( nVertex, ras );
}

bool LayerSurface::GetTargetAtVertex( int nVertex, double* ras )
{
  if ( m_surfaceSource == NULL )
  {
    return false;
  }

  bool bRet = m_surfaceSource->GetRASAtVertex( nVertex, ras );
  if ( bRet )
  {
    //m_volumeRef->RASToTarget( ras, ras );
    m_surfaceSource->ConvertRASToTarget(ras, ras);
  }

  double* offset = GetProperty()->GetPosition();
  for ( int i = 0; i < 3; i++ )
  {
    ras[i] += offset[i];
  }

  return bRet;
}

void LayerSurface::SetActiveSurface( int nSurface )
{
  static int old_thickness = 0;
  if ( m_surfaceSource && m_surfaceSource->SetActiveSurface( nSurface ) )
  {
    if ( GetActiveVector() >= 0 )
    {
      UpdateVectorActor2D();
    }

    double pos[3] = { 0, 0, 0 };
    if (IsInflated() || nSurface == FSSurface::SurfaceInflated)
    {
      pos[0] = (GetHemisphere() == 0 ? -45 : 45);
      GetProperty()->SetPosition(pos);
      int nThickness = GetProperty()->GetEdgeThickness();
      if (nThickness > 0)
        old_thickness = nThickness;
      GetProperty()->SetEdgeThickness(0);
    }
    else
    {
      GetProperty()->SetPosition(pos);
      if (old_thickness > 0)
        GetProperty()->SetEdgeThickness(old_thickness);
    }

    emit ActorUpdated();
  }
}

int LayerSurface::GetActiveSurface()
{
  return ( m_surfaceSource ? m_surfaceSource->GetActiveSurface() : -1 );
}

void LayerSurface::SetActiveVector( int nVector )
{
  if ( m_surfaceSource && m_surfaceSource->SetActiveVector( nVector ) )
  {
    UpdateVectorActor2D();
    SetVisible( IsVisible() );  // refresh visibility
  }
}

int LayerSurface::GetActiveVector()
{
  return ( m_surfaceSource ? m_surfaceSource->GetActiveVector() : -1 );
}

int LayerSurface::GetNumberOfVectorSets()
{
  return ( m_surfaceSource ? m_surfaceSource->GetNumberOfVectorSets() : 0 );
}

void LayerSurface::GetVectorAtVertex( int nVertex, double* vec_out )
{
  if ( m_surfaceSource )
  {
    m_surfaceSource->GetVectorAtVertex( nVertex, vec_out );
  }
}

void LayerSurface::GetNormalAtVertex( int nVertex, double* vec_out )
{
  if ( m_surfaceSource )
  {
    m_surfaceSource->GetNormalAtVertex( nVertex, vec_out );
  }
}

void LayerSurface::GetCurvatureRange( double* range )
{
  if ( m_surfaceSource && m_surfaceSource->GetMRIS() )
  {
    range[0] = m_surfaceSource->GetMRIS()->min_curv;
    range[1] = m_surfaceSource->GetMRIS()->max_curv;
  }
}

double LayerSurface::GetCurvatureValue( int nVertex )
{
  return m_surfaceSource->GetCurvatureValue( nVertex );
}

bool LayerSurface::HasCurvature()
{
  return m_surfaceSource->IsCurvatureLoaded();
}

bool LayerSurface::HasOverlay()
{
  return m_overlays.size() > 0;
}

void LayerSurface::SetActiveOverlay( int nOverlay )
{
  if ( nOverlay < (int)m_overlays.size() )
  {
    if ( m_nActiveOverlay < 0 && nOverlay >= 0 )
    {
      this->GetProperty()->blockSignals(true);
      this->GetProperty()->SetCurvatureMap( LayerPropertySurface::CM_Binary );
      this->GetProperty()->blockSignals(false);
    }
    m_nActiveOverlay = nOverlay;
    UpdateOverlay(false);
    emit ActiveOverlayChanged( nOverlay );
    emit ActorUpdated();
    GetProperty()->SetShowOverlay(true);
  }
}

void LayerSurface::SetActiveOverlay( const QString& name )
{
  for ( int i = 0; i < m_overlays.size(); i++ )
  {
    if ( m_overlays[i]->GetName() == name )
    {
      SetActiveOverlay( i );
      return;
    }
  }
}

int LayerSurface::GetNumberOfOverlays()
{
  return m_overlays.size();
}

SurfaceOverlay* LayerSurface::GetOverlay( const QString& name )
{
  for ( int i = 0; i < m_overlays.size(); i++ )
  {
    if ( m_overlays[i]->GetName() == name )
    {
      return m_overlays[i];
    }
  }

  return NULL;
}

int LayerSurface::GetActiveOverlayIndex()
{
  return m_nActiveOverlay;
}

SurfaceOverlay* LayerSurface::GetActiveOverlay()
{
  if ( m_nActiveOverlay >= 0 )
  {
    return m_overlays[ m_nActiveOverlay ];
  }
  else
  {
    return NULL;
  }
}

SurfaceOverlay* LayerSurface::GetCorrelationOverlay()
{
  for (int i = 0; i < m_overlays.size(); i++)
  {
    if ( m_overlays[i]->HasCorrelationData())
      return m_overlays[i];
  }
  return NULL;
}

SurfaceOverlay* LayerSurface::GetOverlay( int n )
{
  if ( n >= 0 && n < (int)m_overlays.size() )
  {
    return m_overlays[n];
  }
  else
  {
    return NULL;
  }
}

/*
void LayerSurface::CopyCorrelationOverlay(SurfaceOverlay* overlay)
{

}
*/

void LayerSurface::UpdateCorrelationOverlayAtVertex( int nVertex )
{
  SurfaceOverlay* overlay = GetOverlay( m_nActiveOverlay );
  if ( nVertex < 0 || !overlay || !overlay->HasCorrelationData() )
  {
    return;
  }

  overlay->UpdateCorrelationAtVertex( nVertex );
}

void LayerSurface::UpdateCorrelationOverlay()
{
  int nVertex = GetVertexIndexAtTarget( GetSlicePosition(), NULL );
  UpdateCorrelationOverlayAtVertex( nVertex );
}

void LayerSurface::UpdateOverlay( bool bAskRedraw )
{
  vtkPolyDataMapper* mapper = vtkPolyDataMapper::SafeDownCast( m_mainActor->GetMapper() );
  vtkPolyDataMapper* wf_mapper = vtkPolyDataMapper::SafeDownCast( m_wireframeActor->GetMapper() );
  if (!mapper->GetInput() || !wf_mapper->GetInput())
  {
    mapper->Update();
    wf_mapper->Update();
  }
  vtkPolyData* polydata = mapper->GetInput();
  vtkPolyData* polydataWireframe = wf_mapper->GetInput();
  if ( (GetProperty()->GetShowOverlay() && m_nActiveOverlay >= 0) ||
       (GetProperty()->GetShowAnnotation() && m_nActiveAnnotation >= 0) )
  {
    if ( mapper )
    {
      int nCount = polydata->GetPoints()->GetNumberOfPoints();
      vtkSmartPointer<vtkUnsignedCharArray> array = vtkUnsignedCharArray::SafeDownCast( polydata->GetPointData()->GetArray( "Overlay" ) );
      if ( array.GetPointer() == NULL )
      {
        array = vtkSmartPointer<vtkUnsignedCharArray>::New();
        array->SetNumberOfComponents( 4 );
        array->SetNumberOfTuples( nCount );
        array->SetName( "Overlay" );
        polydata->GetPointData()->AddArray( array );
        polydataWireframe->GetPointData()->AddArray( array );
      }
      unsigned char* data = new unsigned char[ nCount*4 ];
      if (polydata->GetPointData()->GetScalars("Curvature"))
        GetProperty()->GetCurvatureLUT()->MapScalarsThroughTable( polydata->GetPointData()->GetScalars("Curvature"), data, VTK_RGBA );
      else
      {
        double* c = GetProperty()->GetBinaryColor();
        int r = (int)(c[0]*255);
        int g = (int)(c[1]*255);
        int b = (int)(c[2]*255);
        for (int i = 0; i < nCount; i++)
        {
          data[i*4] = r;
          data[i*4+1] = g;
          data[i*4+2] = b;
          data[i*4+3] = 255;
        }
      }
      if (GetProperty()->GetShowOverlay() && m_nActiveOverlay >= 0)
        GetActiveOverlay()->MapOverlay( data );
      if (GetProperty()->GetShowAnnotation() && m_nActiveAnnotation >= 0)
        GetActiveAnnotation()->MapAnnotationColor(data);
      MapLabels( data, nCount );
      for ( int i = 0; i < nCount; i++ )
      {
        array->SetTupleValue( i, data + i*4 );
      }
      delete[] data;
      polydata->GetPointData()->SetActiveScalars( "Overlay" );
      if ( GetProperty()->GetMeshColorMap() == LayerPropertySurface::MC_Surface )
      {
        polydataWireframe->GetPointData()->SetActiveScalars( "Overlay" );
      }
    }
  }
  else
  {
    if ( m_labels.size() == 0 || m_nActiveLabel < 0)   // no labels
    {
      polydata->GetPointData()->SetActiveScalars( "Curvature" );
      if ( GetProperty()->GetMeshColorMap() == LayerPropertySurface::MC_Surface )
      {
        polydataWireframe->GetPointData()->SetActiveScalars( "Curvature" );
      }
    }
    else
    {
      if ( mapper )
      {
        int nCount = polydata->GetPoints()->GetNumberOfPoints();
        vtkSmartPointer<vtkUnsignedCharArray> array = vtkUnsignedCharArray::SafeDownCast( polydata->GetPointData()->GetArray( "Overlay" ) );
        if ( array.GetPointer() == NULL )
        {
          array = vtkSmartPointer<vtkUnsignedCharArray>::New();
          array->SetNumberOfComponents( 4 );
          array->SetNumberOfTuples( nCount );
          array->SetName( "Overlay" );
          polydata->GetPointData()->AddArray( array );
          polydataWireframe->GetPointData()->AddArray( array );
        }
        unsigned char* data = new unsigned char[ nCount*4 ];
        if (polydata->GetPointData()->GetScalars("Curvature"))
          GetProperty()->GetCurvatureLUT()->MapScalarsThroughTable( polydata->GetPointData()->GetScalars("Curvature"), data, VTK_RGBA );
        else
        {
          double* dColor = GetProperty()->GetBinaryColor();
          unsigned char rgba[4] = { (int)(dColor[0]*255), (int)(dColor[1]*255), (int)(dColor[2]*255), 255 };
          for (int i = 0; i < nCount*4; i+=4)
            memcpy(data+i, rgba, 4);
        }

        MapLabels( data, nCount );
        for ( int i = 0; i < nCount; i++ )
        {
          array->SetTupleValue( i, data + i*4 );
        }
        delete[] data;
        polydata->GetPointData()->SetActiveScalars( "Overlay" );
        if ( GetProperty()->GetMeshColorMap() == LayerPropertySurface::MC_Surface )
        {
          polydataWireframe->GetPointData()->SetActiveScalars( "Overlay" );
        }
      }
    }
  }
  if ( bAskRedraw )
  {
    emit ActorUpdated();
  }
}

void LayerSurface::UpdateRenderMode()
{
//  m_mainActor->GetProperty()->EdgeVisibilityOff();
//  m_mainActor->GetProperty()->BackfaceCullingOn();
  switch ( GetProperty()->GetSurfaceRenderMode() )
  {
  case LayerPropertySurface::SM_Surface:
    m_mainActor->VisibilityOn();
    m_wireframeActor->VisibilityOff();
    break;
  case LayerPropertySurface::SM_Wireframe:
    m_mainActor->VisibilityOff();
    m_wireframeActor->VisibilityOn();
    m_wireframeActor->GetProperty()->SetLineWidth( 1 );
    break;
  case LayerPropertySurface::SM_SurfaceAndWireframe:
    m_mainActor->VisibilityOn();
    m_wireframeActor->VisibilityOn();
    m_wireframeActor->GetProperty()->SetLineWidth( 2 );
    break;
  }
  emit ActorUpdated();
}

void LayerSurface::SetActiveAnnotation( int n )
{
  if ( n < (int)m_annotations.size() )
  {
    if ( m_nActiveAnnotation < 0 && n >= 0 )
    {
      this->GetProperty()->SetCurvatureMap( LayerPropertySurface::CM_Binary );
    }
    m_nActiveAnnotation = n;
    UpdateOverlay();
    emit ActiveAnnotationChanged( n );
    emit ActorUpdated();
    GetProperty()->SetShowAnnotation(true);
  }
}

void LayerSurface::SetActiveAnnotation( const QString& name )
{
  for ( int i = 0; i < m_annotations.size(); i++ )
  {
    if ( m_annotations[i]->GetName() == name )
    {
      SetActiveAnnotation( i );
      return;
    }
  }
}

int LayerSurface::GetNumberOfAnnotations()
{
  return m_annotations.size();
}

SurfaceAnnotation* LayerSurface::GetAnnotation( const QString& name )
{
  for ( int i = 0; i < m_annotations.size(); i++ )
  {
    if ( m_annotations[i]->GetName() == name )
    {
      return m_annotations[i];
    }
  }

  return NULL;
}

int LayerSurface::GetActiveAnnotationIndex()
{
  return m_nActiveAnnotation;
}

SurfaceAnnotation* LayerSurface::GetActiveAnnotation()
{
  if ( m_nActiveAnnotation >= 0 )
  {
    return m_annotations[ m_nActiveAnnotation ];
  }
  else
  {
    return NULL;
  }
}

SurfaceAnnotation* LayerSurface::GetAnnotation( int n )
{
  if ( n >= 0 && n < (int)m_annotations.size() )
  {
    return m_annotations[n];
  }
  else
  {
    return NULL;
  }
}

void LayerSurface::UpdateVertexRender()
{
  m_vertexActor->SetVisibility( GetProperty()->GetShowVertices()? 1: 0 );
  m_vertexActor->GetProperty()->SetPointSize( GetProperty()->GetVertexPointSize() );
  m_vertexActor->GetProperty()->SetColor( GetProperty()->GetVertexColor() );
  for (int i = 0; i < 3; i++)
  {
    m_vertexActor2D[i]->SetVisibility( GetProperty()->GetShowVertices()? 1: 0 );
    m_vertexActor2D[i]->GetProperty()->SetPointSize( GetProperty()->GetVertexPointSize() );
    m_vertexActor2D[i]->GetProperty()->SetColor( GetProperty()->GetVertexColor() );
  }
  emit ActorUpdated();
}

void LayerSurface::UpdateMeshRender()
{
  vtkPolyDataMapper* mapper = vtkPolyDataMapper::SafeDownCast( m_mainActor->GetMapper() );
  vtkPolyData* polydata = mapper->GetInput();
  vtkPolyDataMapper* mapperWireframe = vtkPolyDataMapper::SafeDownCast( m_wireframeActor->GetMapper() );
  vtkPolyData* polydataWireframe = mapperWireframe->GetInput();
  mapperWireframe->SetScalarVisibility( GetProperty()->GetMeshColorMap() != LayerPropertySurface::MC_Solid ? 1:0 );
  switch ( GetProperty()->GetMeshColorMap() )
  {
  case LayerPropertySurface::MC_Surface:
    if (polydata->GetPointData()->GetScalars())
      polydataWireframe->GetPointData()->SetActiveScalars( polydata->GetPointData()->GetScalars()->GetName() );
    mapperWireframe->SetLookupTable( mapper->GetLookupTable() );
    break;
  case LayerPropertySurface::MC_Curvature:
  {
    // always display as threshold for curvature
    vtkSmartPointer<vtkRGBAColorTransferFunction> lut = vtkSmartPointer<vtkRGBAColorTransferFunction>::New();
    GetProperty()->BuildCurvatureLUT( lut, LayerPropertySurface::CM_Threshold );
    m_wireframeActor->GetMapper()->SetLookupTable( lut );
    polydataWireframe->GetPointData()->SetActiveScalars( "Curvature" );
  }
  break;
  case LayerPropertySurface::MC_Overlay:
    polydataWireframe->GetPointData()->SetActiveScalars( "Overlay" );
    break;
  case LayerPropertySurface::MC_Solid:
    m_wireframeActor->GetProperty()->SetColor( GetProperty()->GetMeshColor() );
    break;
  default:
    break;
  }
  emit ActorUpdated();
}

void LayerSurface::UpdateActorPositions()
{
  double* pos = GetProperty()->GetPosition();
  for ( int i = 0; i < 3; i++ )
  {
    m_sliceActor3D[i]->SetPosition( pos );
  }

  m_sliceActor2D[0]->SetPosition( 0.1, pos[1], pos[2] );
  m_sliceActor2D[1]->SetPosition( pos[0], 0.1, pos[2] );
  m_sliceActor2D[2]->SetPosition( pos[0], pos[1], -0.1 );
  m_vectorActor2D[0]->SetPosition( 1.0, pos[1], pos[2] );
  m_vectorActor2D[1]->SetPosition( pos[0], 1.0, pos[2] );
  m_vectorActor2D[2]->SetPosition( pos[0], pos[1], -1.0 );

  m_vertexActor2D[0]->SetPosition( 0.1, pos[1], pos[2] );
  m_vertexActor2D[1]->SetPosition( pos[0], 0.1, pos[2] );
  m_vertexActor2D[2]->SetPosition( pos[0], pos[1], -0.1 );

  m_mainActor->SetPosition( pos );
  m_vectorActor->SetPosition( pos );
  m_vertexActor->SetPosition( pos );
  m_wireframeActor->SetPosition( pos );

  emit ActorUpdated();
}

void LayerSurface::UpdateROIPosition(double dx, double dy, double dz)
{
  m_roi->GetActor()->AddPosition(dx, dy, dz);
}

int LayerSurface::GetHemisphere()
{
  return m_surfaceSource->GetMRIS()->hemisphere;
}

int LayerSurface::GetNumberOfLabels()
{
  return m_labels.size();
}

SurfaceLabel* LayerSurface::GetLabel( int n )
{
  if ( n >= 0 && n < (int)m_labels.size() )
  {
    return m_labels[n];
  }
  else
  {
    return NULL;
  }
}

void LayerSurface::SetActiveLabel( int n )
{
  if ( n < (int)m_labels.size() && n != m_nActiveLabel )
  {
    m_nActiveLabel = n;
    /*
    UpdateColorMap();
    for (int i = 0; i < m_labels.size(); i++)
    {
      m_labels[i]->GetOutlineActor()->VisibilityOff();
    }
    if (n >= 0 && m_labels[n]->GetShowOutline())
      m_labels[n]->GetOutlineActor()->VisibilityOn();

    emit ActorUpdated();
    */
    emit ActiveLabelChanged( n );
  }
}

void LayerSurface::MoveLabelToTop(SurfaceLabel *label)
{
  SurfaceLabel* activeLabel = GetActiveLabel();
  for (int i = 0; i < m_labels.size(); i++)
  {
    if (label == m_labels[i])
    {
      m_labels.removeAt(i);
      m_labels.insert(0, label);
      SetActiveLabel(activeLabel);

      UpdateOverlay(false);
      emit Modified();
      emit ActorChanged();
      return;
    }
  }
}

void LayerSurface::DeleteLabel(SurfaceLabel *label)
{
  SurfaceLabel* activeLabel = GetActiveLabel();
  for (int i = 0; i < m_labels.size(); i++)
  {
    if (label == m_labels[i])
    {
      m_labels.removeAt(i);
      if (activeLabel == label)
      {
        if (m_labels.isEmpty())
          m_nActiveLabel = -1;
        else
        {
          if (i >= m_labels.size())
            SetActiveLabel(m_labels.size()-1);
          else
            SetActiveLabel(i);
        }
      }
      delete label;
      emit SurfaceLabelDeleted(label);

      UpdateOverlay(false);
      emit Modified();
      emit ActorChanged();
      return;
    }
  }
}

void LayerSurface::SetActiveLabel(SurfaceLabel *label)
{
  for (int i = 0; i < m_labels.size(); i++)
  {
    if (label == m_labels[i])
    {
      SetActiveLabel(i);
      return;
    }
  }
}

void LayerSurface::MapLabels( unsigned char* data, int nVertexCount )
{
  for ( int i = m_labels.size()-1; i >= 0; i-- )
  {
    if (m_labels[i]->IsVisible())
      m_labels[i]->MapLabel( data, nVertexCount );
  }
}

void LayerSurface::SetActiveLabelColor(const QColor &c)
{
  if ( m_nActiveLabel >= 0)
  {
    m_labels[m_nActiveLabel]->SetColor(c.redF(), c.greenF(), c.blueF());
    UpdateColorMap();
    emit ActorUpdated();
  }
}

void LayerSurface::SetActiveLabelOutline(bool bOutline)
{
  if ( m_nActiveLabel >= 0)
  {
    m_labels[m_nActiveLabel]->SetShowOutline(bOutline);

    UpdateColorMap();
    emit ActorUpdated();
  }
}

void LayerSurface::SetActiveAnnotationOutline(bool bOutline)
{
  if ( m_nActiveAnnotation >= 0)
  {
    m_annotations[m_nActiveAnnotation]->SetShowOutline(bOutline);
    UpdateColorMap();
    emit ActorUpdated();
  }
}

void LayerSurface::RepositionSurface( LayerMRI* mri, int nVertex, double value, int size, double sigma, int flags )
{
  m_surfaceSource->Reposition( mri->GetSourceVolume(), nVertex, value, size, sigma, flags );
  SetModified();
  m_bUndoable = true;
  emit ActorUpdated();
}

void LayerSurface::RepositionSurface( LayerMRI* mri, int nVertex, double* pos, int size, double sigma, int flags )
{
  m_surfaceSource->Reposition( mri->GetSourceVolume(), nVertex, pos, size, sigma, flags );
  SetModified();
  m_bUndoable = true;
  emit ActorUpdated();
}

void LayerSurface::RepositionSmoothSurface(int nVertex, int size, int n_steps)
{
  m_surfaceSource->RepositionSmooth(nVertex, size, n_steps);
  SetModified();
  m_bUndoable = true;
  emit ActorUpdated();
}

bool LayerSurface::SmoothSurface(int nMethod, int niters, double lambda, double k_cutoff)
{
  bool ret = m_surfaceSource->Smooth(nMethod, niters, lambda, k_cutoff);
  if (ret)
  {
    SetModified();
    m_bUndoable = true;
    emit ActorUpdated();
  }
  return ret;
}

void LayerSurface::RemoveIntersections()
{
  m_surfaceSource->RemoveIntersections();
  SetModified();
  m_bUndoable = true;
  emit ActorUpdated();
}

void LayerSurface::RepositionVertex(int nVertex, double* pos)
{
  m_surfaceSource->RepositionVertex(nVertex, pos);
  SetModified();
  m_bUndoable = true;
  emit ActorUpdated();
}

void LayerSurface::Undo()
{
  m_surfaceSource->UndoReposition();
  SetModified();
  m_bUndoable = false;
  emit ActorUpdated();
}

bool LayerSurface::HasUndo()
{
  return m_bUndoable;
}

bool LayerSurface::HasValidVolumeGeometry()
{
  return m_surfaceSource && m_surfaceSource->HasValidVolumeGeometry();
}

int LayerSurface::GetNumberOfVertices()
{
  return this->m_surfaceSource->GetNumberOfVertices();
}

void LayerSurface::ResetVolumeRef()
{
  m_volumeRef = NULL;
  m_surfaceSource->ResetVolumeRef();
}

void LayerSurface::SetCurrentVertex(int n)
{
  if (m_nCurrentVertex != n)
  {
    m_nCurrentVertex = n;
    emit CurrentVertexChanged(n);
  }
}

bool LayerSurface::GetCorrelationOverlayDataAtVertex(int nVert, float *output, int nFrames)
{
  SurfaceOverlay* overlay = NULL;
  for (int i = 0; i < m_overlays.size(); i++)
  {
    if (m_overlays[i]->GetNumberOfFrames() == nFrames)
    {
      overlay = m_overlays[i];
      break;
    }
  }

  if (overlay)
  {
    return overlay->GetDataAtVertex(nVert, output);
  }
  return false;
}

bool LayerSurface::IsInflated()
{
  return GetFileName().toLower().contains("inflated");
}

bool LayerSurface::GetActiveLabelCentroidPosition(double *pos)
{
  SurfaceLabel* label = GetActiveLabel();
  int nvo;
  double x, y, z;
  if (label && label->GetCentroid(&x, &y, &z, &nvo))
  {
    return GetTargetAtVertex(nvo, pos);
  }
}
