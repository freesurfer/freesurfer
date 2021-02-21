/**
 * @brief Layer data object for MRI surface.
 *
 */
/*
 * Original Author: Ruopeng Wang
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
#include "vtkCellArray.h"
#include "vtkSTLWriter.h"
#include "vtkTransformPolyDataFilter.h"
#include "vtkTransform.h"
#include "FSSurface.h"
#include "LayerMRI.h"
#include "SurfaceAnnotation.h"
#include "SurfaceLabel.h"
#include "LayerPropertySurface.h"
#include "SurfaceROI.h"
#include <QFileInfo>
#include <QDir>
#include <QTextStream>
#include <QFile>
#include <QDebug>
#include "SurfaceOverlayProperty.h"
#include "MyUtils.h"
#include "SurfaceOverlay.h"
#include "SurfaceSpline.h"
#include "vtkMaskPoints.h"
#include "vtkExtractPolyDataGeometry.h"
#include "vtkBox.h"
#include "vtkDoubleArray.h"
#include "LayerROI.h"
#include "LayerPropertyROI.h"
#include "SurfacePath.h"
#include <QSet>

LayerSurface::LayerSurface( LayerMRI* ref, QObject* parent ) : LayerEditable( parent ),
  m_surfaceSource( NULL ),
  m_bResampleToRAS( true ),
  m_volumeRef( ref ),
  m_nActiveOverlay( -1 ),
  m_nActiveAnnotation( -1 ),
  m_nActiveLabel( -1 ),
  m_nActiveSpline(-1),
  m_bUndoable( false ),
  m_bVector2DPendingUpdate( true ),
  m_bLoadAll(false),
  m_nCurrentVertex(-1),
  m_bVisibleIn3D(true),
  m_nActiveRGBMap(-1),
  m_nColorDataCache(NULL),
  m_surfaceContralateral(NULL),
  m_surfaceSphere1(NULL),
  m_surfaceSphere2(NULL),
  m_nMouseVertex(-1),
  m_nActivePath(-1),
  m_marks(NULL)
{
  m_strTypeNames.push_back( "Surface" );
  m_sPrimaryType = "Surface";

  m_sMappingSurfaceName = "white";
  m_sSphereFilename = "sphere.reg";

  // create property before actors!
  mProperty = new LayerPropertySurface( this );

  for ( int i = 0; i < 3; i++ )
  {
    // m_nSliceNumber[i] = 0;
    m_sliceActor2D[i] = vtkSmartPointer<vtkActor>::New();
    m_sliceActor3D[i] = vtkSmartPointer<vtkActor>::New();
    m_vectorActor2D[i] = vtkSmartPointer<vtkActor>::New();
    m_vertexActor2D[i] = vtkSmartPointer<vtkActor>::New();
    m_vertexActor2D[i]->SetProperty( m_vertexActor2D[i]->MakeProperty() );
    m_vertexActor2D[i]->GetProperty()->SetRepresentationToPoints();
    m_vertexActor2D[i]->VisibilityOff();
  }

  m_mainActor = vtkSmartPointer<vtkActor>::New();
  m_mainActor->GetProperty()->SetEdgeColor( 0.75, 0.75, 0.75 );

  double ratio = 1;
#if VTK_MAJOR_VERSION > 7
  ratio = MainWindow::GetMainWindow()->devicePixelRatio();
#endif
  m_vectorActor = vtkSmartPointer<vtkActor>::New();
  m_vectorActor->GetProperty()->SetColor( GetProperty()->GetVectorColor() );
  m_vectorActor->GetProperty()->SetPointSize( GetProperty()->GetVectorPointSize()*ratio );
  m_vectorActor->PickableOff();

  for ( int i = 0; i < 3; i++ )
  {
    m_vectorActor2D[i]->GetProperty()->SetColor( GetProperty()->GetVectorColor() );
    m_vectorActor2D[i]->GetProperty()->SetPointSize( GetProperty()->GetVectorPointSize()*ratio );
    m_vectorActor2D[i]->PickableOff();
  }

  m_vertexActor = vtkSmartPointer<vtkActor>::New();
  m_vertexActor->GetProperty()->SetRepresentationToPoints();
  m_vertexActor->VisibilityOff();

  m_wireframeActor = vtkSmartPointer<vtkActor>::New();
  m_wireframeActor->VisibilityOff();

  m_roi = new SurfaceROI(this);

  //  m_spline = new SurfaceSpline(this);
  //  connect(m_spline, SIGNAL(SplineChanged()), this, SIGNAL(ActorChanged()));

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

#if VTK_MAJOR_VERSION > 5
  m_mainActor->ForceOpaqueOn();
  m_wireframeActor->ForceOpaqueOn();
  m_vectorActor->ForceOpaqueOn();
  for (int i = 0; i < 3; i++)
    m_vectorActor2D[i]->ForceOpaqueOn();
#endif
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

  for (int i = 0; i < m_paths.size(); i++)
    delete m_paths[i];

  if (m_marks)
    delete m_marks;

  if (m_nColorDataCache)
    delete[] m_nColorDataCache;
}

void LayerSurface::SetRefVolume(LayerMRI *ref)
{
  m_volumeRef = ref;
  if (m_volumeRef)
    connect( m_volumeRef, SIGNAL(destroyed()), this, SLOT(ResetVolumeRef()), Qt::UniqueConnection);
}

bool LayerSurface::LoadSurfaceFromFile(bool bIgnoreVG)
{
  if ( m_surfaceSource )
  {
    delete m_surfaceSource;
  }

  m_surfaceSource = new FSSurface( m_volumeRef ? m_volumeRef->GetSourceVolume() : NULL );
  m_surfaceSource->SetIgnoreVolumeGeometry(bIgnoreVG);
  setProperty("IgnoreVG", bIgnoreVG);
  if ( !m_surfaceSource->MRISRead( m_sFilename,
                                   m_sVectorFilename,
                                   m_sPatchFilename,
                                   m_sTargetFilename,
                                   m_sSphereFilename,
                                   m_listSupFiles)
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
  m_nColorDataCache = new unsigned char[GetNumberOfVertices()*4];
}

bool LayerSurface::CreateFromMRIS(void *mris_ptr)
{
  MRIS* mris = (MRIS*)mris_ptr;
  m_surfaceSource = new FSSurface( m_volumeRef ? m_volumeRef->GetSourceVolume() : NULL );
  if ( !m_surfaceSource->CreateFromMRIS(mris) )
  {
    return false;
  }
  SetFileName(mris->fname.data());
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

bool LayerSurface::SaveSurfaceAsSTL(const QString &fn)
{
  vtkSmartPointer<vtkTransform> tr = vtkSmartPointer<vtkTransform>::New();
  tr->DeepCopy(m_surfaceSource->GetSurfaceToRasTransform());
  vtkSmartPointer<vtkTransformPolyDataFilter> filter =
        vtkSmartPointer<vtkTransformPolyDataFilter>::New();
  filter->SetTransform( tr );
#if VTK_MAJOR_VERSION > 5
  filter->SetInputData( m_surfaceSource->GetPolyData() );
#else
  filter->SetInput( m_surfaceSource->GetPolyData() );
#endif
  filter->Update();
  vtkSmartPointer<vtkSTLWriter> writer = vtkSmartPointer<vtkSTLWriter>::New();
  writer->SetFileName(fn.toLatin1().constData());
#if VTK_MAJOR_VERSION > 5
  writer->SetInputData( filter->GetOutput() );
#else
  writer->SetInput( filter->GetOutput() );
#endif
  return writer->Write();
}

bool LayerSurface::SaveSurface( )
{
  if ( m_sFilename.size() == 0 )
  {
    cerr << "No filename provided to save surface.\n";
    return false;
  }

  return SaveSurface( m_sFilename.toLatin1().data() );
}

bool LayerSurface::WriteIntersection(const QString &filename, int nPlane, LayerMRI* mri_ref)
{
  vtkPolyData* polydata = vtkPolyData::SafeDownCast(m_sliceActor2D[nPlane]->GetMapper()->GetInput());
  if (!polydata)
  {
    return false;
  }

  vtkCellArray* lines = polydata->GetLines();
  vtkIdType nPts, *pts;
  vtkPoints* points = polydata->GetPoints();

  if (lines && points)
  {
    QFile file(filename);
    if (!file.open(QIODevice::WriteOnly | QIODevice::Text))
    {
      cerr << "Cannot write file " << qPrintable(filename) << endl;
      return false;
    }
    QTextStream out(&file);
    QString str[3] = {"sag", "cor", "hor"};
    out << "#orientation " << str[nPlane] << "\n";
    out << "#points " << points->GetNumberOfPoints() << "\n";
    for (int i = 0; i < points->GetNumberOfPoints(); i++)
    {
      double* pt = points->GetPoint(i);
      double ras[3];
      mri_ref->TargetToRAS(pt, ras);
      mri_ref->RASToOriginalIndex(ras, ras);
      out << ras[0] << " " << ras[1] << " " << ras[2] << "\n";
    }
    out << "#lines " << lines->GetNumberOfCells() << "\n";
    lines->InitTraversal();
    while (lines->GetNextCell(nPts, pts))
    {
      for (int i = 0; i < nPts; i++)
        out << pts[i] << " ";
      out << "\n";
    }
    cout << "Intersection data written to " << qPrintable(filename) << "\n";
  }
  return true;
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

  GetProperty()->RebuildCurvatureLUT();
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
  overlay->SetRegFileName( fn_reg );

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
  if (overlay->HasCorrelationData())
    UpdateCorrelationOverlay();
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
  connect(annot, SIGNAL(Modified()), SLOT(UpdateColorMap()));

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

  InitializeLabel(label);
  return true;
}

void LayerSurface::InitializeLabel(SurfaceLabel *label)
{
  m_labels.insert(0, label);
  connect(label, SIGNAL(SurfaceLabelChanged()), this, SLOT(UpdateColorMap()));
  connect(label, SIGNAL(SurfaceLabelChanged()), GetProperty(), SIGNAL(PropertyChanged()));
  connect(label, SIGNAL(SurfaceLabelVisibilityChanged()), this, SLOT(UpdateColorMap()));

  SetActiveLabel( 0 );

  UpdateOverlay(false);

  emit Modified();
  emit SurfaceLabelAdded( label );
  emit ActorChanged();
}

SurfaceLabel* LayerSurface::CreateNewLabel(const QString& name_in)
{
  static int ncount = 1;
  QString name = name_in;
  if (name.isEmpty())
  {
    name = QString("label_%1").arg(ncount++);
  }
  SurfaceLabel* label = new SurfaceLabel( this, true );
  InitializeLabel(label);
  label->SetName(name);
  return label;
}

bool LayerSurface::LoadSplineFromFile(const QString &filename)
{
  SurfaceSpline* spline = new SurfaceSpline(this);
  if (!spline->Load(filename))
  {
    spline->deleteLater();
    return false;
  }
  m_splines.insert(0, spline);
  connect(spline, SIGNAL(SplineChanged()), this, SLOT(UpdateColorMap()));
  SetActiveSpline(0);

  emit ActorChanged();
  emit SurfaceSplineAdded( spline );
  return true;
}

void LayerSurface::DeleteSpline(SurfaceSpline *spline)
{
  SurfaceSpline* activeSpline = GetActiveSpline();
  for (int i = 0; i < m_splines.size(); i++)
  {
    if (spline == m_splines[i])
    {
      m_splines.removeAt(i);
      if (activeSpline == spline)
      {
        if (m_splines.isEmpty())
          m_nActiveSpline = -1;
        else
        {
          if (i >= m_splines.size())
            SetActiveSpline(m_splines.size()-1);
          else
            SetActiveSpline(i);
        }
      }
      spline->deleteLater();
      emit SurfaceSplineDeleted(spline);
      emit Modified();
      emit ActorChanged();
      return;
    }
  }
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
#if VTK_MAJOR_VERSION > 5
  mapper->SetInputData( m_surfaceSource->GetPolyData() );
#else
  mapper->SetInput( m_surfaceSource->GetPolyData() );
#endif
  m_mainActor->SetMapper( mapper );
  mapper->Update();

  // vector actor
  mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  vtkSmartPointer<vtkTubeFilter> tube = vtkSmartPointer<vtkTubeFilter>::New();
#if VTK_MAJOR_VERSION > 5
  tube->SetInputData(m_surfaceSource->GetVectorPolyData());
#else
  tube->SetInput(m_surfaceSource->GetVectorPolyData());
#endif
  tube->SetNumberOfSides(8);
  tube->SetRadius(0.04);
  tube->CappingOn();
  mapper->SetInputConnection( tube->GetOutputPort() );
  m_vectorActor->SetMapper( mapper );
  //  mapper->Update();

  for ( int i = 0; i < 3; i++ )
  {
    mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
#if VTK_MAJOR_VERSION > 5
    mapper->SetInputData(  m_surfaceSource->GetVector2DPolyData( i ) );
#else
    mapper->SetInput(  m_surfaceSource->GetVector2DPolyData( i ) );
#endif
    m_vectorActor2D[i]->SetMapper( mapper );
  }

  // vertex actor
  mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
#if VTK_MAJOR_VERSION > 5
  mapper->SetInputData(  m_surfaceSource->GetVertexPolyData() );
#else
  mapper->SetInput(  m_surfaceSource->GetVertexPolyData() );
#endif
  m_vertexActor->SetMapper( mapper );

  // wireframe actor
  mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
#if VTK_MAJOR_VERSION > 5
  mapper->SetInputData(  m_surfaceSource->GetWireframePolyData() );
#else
  mapper->SetInput(  m_surfaceSource->GetWireframePolyData() );
#endif
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
#if VTK_MAJOR_VERSION > 5
    extract->SetInputData(m_surfaceSource->GetPolyData());
#else
    extract->SetInput(m_surfaceSource->GetPolyData());
#endif
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
    vtkSmartPointer<vtkPolyDataMapper> mapper2 = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper2->SetInputConnection( m_cutter[i]->GetOutputPort() );
    //
    // Actors in the scene, drawing the mapped lines.
    //

    double line_w = GetProperty()->GetEdgeThickness();
    double ratio = 1;
#if VTK_MAJOR_VERSION > 7
    ratio = MainWindow::GetMainWindow()->devicePixelRatio();
#endif
    m_sliceActor2D[i]->SetMapper( mapper );
    //  m_sliceActor2D[i]->SetBackfaceProperty( m_sliceActor2D[i]->MakeProperty() );
    //  m_sliceActor2D[i]->GetBackfaceProperty()->BackfaceCullingOff();
    m_sliceActor2D[i]->SetProperty( m_sliceActor2D[i]->MakeProperty() );
    m_sliceActor2D[i]->GetProperty()->SetInterpolationToFlat();
    m_sliceActor2D[i]->GetProperty()->SetLineWidth( line_w*ratio );

    m_sliceActor3D[i]->SetMapper( mapper2 );
    //  m_sliceActor3D[i]->SetBackfaceProperty( m_sliceActor3D[i]->MakeProperty() );
    //  m_sliceActor3D[i]->GetBackfaceProperty()->BackfaceCullingOff();
    m_sliceActor3D[i]->SetProperty( m_sliceActor3D[i]->MakeProperty() );
    m_sliceActor3D[i]->GetProperty()->SetLineWidth( line_w*ratio );
    //    m_sliceActor3D[i]->GetProperty()->SetInterpolationToFlat();

    vtkSmartPointer<vtkPolyDataMapper> mapper3 = vtkSmartPointer<vtkPolyDataMapper>::New();
    m_vertexPoly2D[i] = vtkSmartPointer<vtkPolyData>::New();
#if VTK_MAJOR_VERSION > 5
    mapper3->SetInputData(m_vertexPoly2D[i]);
#else
mapper3->SetInput(m_vertexPoly2D[i]);
#endif 

    mapper3->ScalarVisibilityOff();
    m_vertexActor2D[i]->SetMapper(mapper3);
    m_vertexActor2D[i]->SetProperty( m_vertexActor2D[i]->MakeProperty() );
    m_vertexActor2D[i]->GetProperty()->SetPointSize(3*ratio);
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
  double line_w = GetProperty()->GetEdgeThickness();
#if VTK_MAJOR_VERSION > 7
  line_w *= MainWindow::GetMainWindow()->devicePixelRatio();
#endif
  for ( int i = 0; i < 3; i++ )
  {
    m_sliceActor2D[i]->GetProperty()->SetLineWidth( line_w );
    m_sliceActor3D[i]->GetProperty()->SetLineWidth( line_w );
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
  double ratio = 1;
#if VTK_MAJOR_VERSION > 7
  ratio = MainWindow::GetMainWindow()->devicePixelRatio();
#endif
  for ( int i = 0; i < 3; i++ )
  {
    m_vectorActor2D[i]->GetProperty()->SetPointSize( GetProperty()->GetVectorPointSize()*ratio );
  }
  m_vectorActor->GetProperty()->SetPointSize( GetProperty()->GetVectorPointSize()*ratio );
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
    m_sliceActor2D[i]->GetMapper()->SetScalarVisibility(GetProperty()->GetUseSurfaceColorOn2D()?1:0);
    m_sliceActor3D[i]->GetMapper()->SetScalarVisibility(GetProperty()->GetUseSurfaceColorOn2D()?1:0);
    m_sliceActor2D[i]->GetMapper()->SetLookupTable( GetProperty()->GetCurvatureLUT() );
    m_sliceActor3D[i]->GetMapper()->SetLookupTable( GetProperty()->GetCurvatureLUT() );
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
  foreach (SurfaceSpline* spline, m_splines)
    spline->AppendProp2D(renderer, nPlane);
}

void LayerSurface::Append3DProps( vtkRenderer* renderer, bool* bSliceVisibility )
{
  if (!m_bVisibleIn3D)
    return;

  for ( int i = 0; i < 3; i++ )
  {
    if (bSliceVisibility[i])
    {
      //  renderer->AddViewProp( m_sliceActor3D[i] );
    }
  }

  renderer->AddViewProp( m_mainActor );
  renderer->AddViewProp( m_vectorActor );
  renderer->AddViewProp( m_vertexActor );
  renderer->AddViewProp( m_wireframeActor );

  foreach (SurfaceSpline* spline, m_splines)
    spline->AppendProp3D(renderer);

  foreach (SurfacePath* path, m_paths)
    path->AppendProps(renderer);

  if (m_marks)
    m_marks->AppendProps(renderer);

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
  double dMinVS = qMin(qMin(m_dWorldVoxelSize[0], m_dWorldVoxelSize[1]), m_dWorldVoxelSize[2]);
  switch ( nPlane )
  {
  case 0:
    mReslicePlane[0]->SetOrigin( m_dSlicePosition[0]-pos[0], 0, 0  );
    m_sliceActor2D[0]->SetPosition( 0.1, pos[1], pos[2] );
    m_vertexActor2D[0]->SetPosition( dMinVS+0.1, pos[1], pos[2] );
    m_vectorActor2D[0]->SetPosition( 1.0, pos[1], pos[2] );
    break;
  case 1:
    mReslicePlane[1]->SetOrigin( 0, m_dSlicePosition[1]-pos[1], 0 );
    m_sliceActor2D[1]->SetPosition( pos[0], 0.1, pos[2] );
    m_vertexActor2D[1]->SetPosition( pos[0], dMinVS+0.1, pos[2] );
    m_vectorActor2D[1]->SetPosition( pos[0], 1.0, pos[2] );
    break;
  case 2:
    mReslicePlane[2]->SetOrigin( 0, 0, m_dSlicePosition[2]-pos[2]  );
    m_sliceActor2D[2]->SetPosition( pos[0], pos[1], -0.1 );
    m_vertexActor2D[2]->SetPosition( pos[0], pos[1], -0.1-dMinVS );
    m_vectorActor2D[2]->SetPosition( pos[0], pos[1], -1.0 );
    break;
  }
  double dLen = m_surfaceSource->GetMaxSegmentLength();
  bounds[nPlane*2] = m_dSlicePosition[nPlane]-dLen/2;
  bounds[nPlane*2+1] = bounds[nPlane*2]+dLen/2;
  m_box[nPlane]->SetBounds(bounds);

  dLen = dMinVS/5;
  if (GetProperty()->GetShowVertices())
  {
      vtkPoints* all_pts = m_surfaceSource->GetPolyData()->GetPoints();
      vtkSmartPointer<vtkPoints> pts = vtkSmartPointer<vtkPoints>::New();
      double pt[3];
      for (vtkIdType i = 0; i < all_pts->GetNumberOfPoints(); i++)
      {
          all_pts->GetPoint(i, pt);
          if (qAbs(m_dSlicePosition[nPlane]-pt[nPlane]) < dLen)
          {
              pts->InsertNextPoint(pt);
          }
      }
      vtkSmartPointer<vtkCellArray> verts = vtkSmartPointer<vtkCellArray>::New();
      verts->Allocate( pts->GetNumberOfPoints() );
      for ( int i = 0; i < pts->GetNumberOfPoints(); i++ )
      {
        vtkIdType n = i;
        verts->InsertNextCell( 1, &n );
      }
      m_vertexPoly2D[nPlane]->SetPoints(pts);
      m_vertexPoly2D[nPlane]->SetVerts(verts);
  }

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
}

void LayerSurface::OnSlicePositionChanged3D()
{
  if (!m_splines.isEmpty() && MainWindow::GetMainWindow()->GetSplinePicking())
  {
    int nVertex = this->GetVertexIndexAtTarget(m_dSlicePosition, NULL);
    foreach (SurfaceSpline* spline, m_splines)
    {
      if (!spline->IsLocked())
        spline->SetActiveVertex(nVertex);
    }
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

  foreach (SurfaceSpline* spline, m_splines)
  {
    if (spline->IsVisible() && bVisible)
      spline->SetActorVisible(true);
    else if (!bVisible)
      spline->SetActorVisible(false);
  }

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

int LayerSurface::GetVertexIndexAtTarget( double* pos, double* distance, int surface_type )
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
  return m_surfaceSource->FindVertexAtRAS( realRas, distance, surface_type );
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

void LayerSurface::GetSurfaceRASAtRAS(double* ras_in, double* tkras_out)
{
  if ( m_surfaceSource == NULL )
  {
    return;
  }
  m_surfaceSource->ConvertRASToSurface( ras_in, tkras_out );
}

void LayerSurface::GetRASAtTarget(double *pos_in, double *ras_out)
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
  m_surfaceSource->ConvertTargetToRAS( pos_o, ras_out );
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

int LayerSurface::GetVertexAtSurfaceRAS(double *ras, double *distance)
{
  if (m_surfaceSource == NULL)
    return -1;

  return m_surfaceSource->FindVertexAtSurfaceRAS(ras, distance);
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

    for (int i = 0; i < m_paths.size(); i++)
    {
      m_paths[i]->Update();
    }

    if (m_marks)
      m_marks->Update();

    emit ActorUpdated();
    emit ActiveSurfaceChanged(nSurface);
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
      if (this->GetProperty()->GetCurvatureMap() == LayerPropertySurface::CM_Threshold)
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
    {        if (this->GetProperty()->GetCurvatureMap() == LayerPropertySurface::CM_Threshold)
        this->GetProperty()->SetCurvatureMap( LayerPropertySurface::CM_Binary );
      SetActiveOverlay( i );
      return;
    }        if (this->GetProperty()->GetCurvatureMap() == LayerPropertySurface::CM_Threshold)
      this->GetProperty()->SetCurvatureMap( LayerPropertySurface::CM_Binary );
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
  int nVertex = GetCurrentVertex();
  if (!IsInflated())
    nVertex = GetVertexIndexAtTarget( GetSlicePosition(), NULL );
  UpdateCorrelationOverlayAtVertex( nVertex );
}

void LayerSurface::UpdateOverlay(bool bAskRedraw, bool pre_cached)
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
       (GetProperty()->GetShowAnnotation() && m_nActiveAnnotation >= 0))
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

      QMultiMap<int, QString> zorders;
      zorders.insert(GetProperty()->GetZOrderAnnotation(), "annot");
      zorders.insert(GetProperty()->GetZOrderLabel(), "label");
      zorders.insert(GetProperty()->GetZOrderOverlay(), "overlay");
      bool bLabelOnTop = (zorders.values(zorders.keys().last()).last() == "label");

      unsigned char* data = new unsigned char[ nCount*4 ];
      if (pre_cached && bLabelOnTop)
      {
        memcpy(data, m_nColorDataCache, nCount*4);
        MapLabels( data, nCount );
      }
      else
      {
        if (polydata->GetPointData()->GetScalars("Curvature") && m_nActiveRGBMap < 0)
        {
          GetProperty()->GetCurvatureLUT()->MapScalarsThroughTable( polydata->GetPointData()->GetScalars("Curvature"), data, VTK_RGBA );
        }
        else
        {
          if (m_nActiveRGBMap < 0)
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
          else
          {
            QList<int>& rgb = m_rgbMaps[m_nActiveRGBMap].data;
            for (int i = 0; i < nCount; i++)
            {
              data[i*4] = rgb[i*3];
              data[i*4+1] = rgb[i*3+1];
              data[i*4+2] = rgb[i*3+2];
              data[i*4+3] = 255;
            }
          }
        }
        QList<int> zkeys = zorders.keys();
        for (int i = 0; i < zkeys.size(); i++)
        {
          QStringList values = zorders.values(zkeys[i]);
          for (int j = 0; j < values.size(); j++)
          {
            QString render = values[j];
            if ( render == "overlay")
            {
              if (GetProperty()->GetShowOverlay() && m_nActiveOverlay >= 0)
                GetActiveOverlay()->MapOverlay( data );
            }
            else if (render == "annot")
            {
              if (GetProperty()->GetShowAnnotation() && m_nActiveAnnotation >= 0)
                GetActiveAnnotation()->MapAnnotationColor(data);
            }
            else if (render == "label")
            {
              if (bLabelOnTop)
                memcpy(m_nColorDataCache, data, nCount*4);
              MapLabels( data, nCount );
            }
          }
        }
      }
      for ( int i = 0; i < nCount; i++ )
      {
#if VTK_MAJOR_VERSION > 5
        array->SetTypedTuple( i, data + i*4 );
#else
        array->SetTupleValue( i, data + i*4 );
#endif
      }
      array->Modified();
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
    if ( ((m_labels.isEmpty() && m_mappedLabels.isEmpty()) || (m_mappedLabels.isEmpty() && m_nActiveLabel < 0)) && m_nActiveRGBMap < 0)   // no labels
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
        vtkIdType nCount = polydata->GetPoints()->GetNumberOfPoints();
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
        if (pre_cached)
          memcpy(data, m_nColorDataCache, nCount*4);
        else
        {
          if (polydata->GetPointData()->GetScalars("Curvature") && m_nActiveRGBMap < 0)
            GetProperty()->GetCurvatureLUT()->MapScalarsThroughTable( polydata->GetPointData()->GetScalars("Curvature"), data, VTK_RGBA );
          else
          {
            if (m_nActiveRGBMap < 0)
            {
              double* dColor = GetProperty()->GetBinaryColor();
              unsigned char rgba[4] = { (unsigned char)(dColor[0]*255), (unsigned char)(dColor[1]*255), (unsigned char)(dColor[2]*255), 255 };
              for (vtkIdType i = 0; i < nCount*4; i+=4)
                memcpy(data+i, rgba, 4);
            }
            else
            {
              QList<int>& rgb = m_rgbMaps[m_nActiveRGBMap].data;
              for (vtkIdType i = 0; i < nCount; i++)
              {
                data[i*4] = rgb[i*3];
                data[i*4+1] = rgb[i*3+1];
                data[i*4+2] = rgb[i*3+2];
                data[i*4+3] = 255;
              }
            }
          }
          memcpy(m_nColorDataCache, data, nCount*4);
        }

        MapLabels( data, nCount );
        for ( int i = 0; i < nCount; i++ )
        {
#if VTK_MAJOR_VERSION > 5
          array->SetTypedTuple( i, data + i*4 );
#else
          array->SetTupleValue( i, data + i*4 );
#endif
        }
        array->Modified();
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
  double line_w = 1;
#if VTK_MAJOR_VERSION > 7
  line_w = MainWindow::GetMainWindow()->devicePixelRatio();
#endif
  switch ( GetProperty()->GetSurfaceRenderMode() )
  {
  case LayerPropertySurface::SM_Surface:
    m_mainActor->VisibilityOn();
    m_wireframeActor->VisibilityOff();
    break;
  case LayerPropertySurface::SM_Wireframe:
    m_mainActor->VisibilityOff();
    m_wireframeActor->VisibilityOn();
    m_wireframeActor->GetProperty()->SetLineWidth( line_w );
    break;
  case LayerPropertySurface::SM_SurfaceAndWireframe:
    m_mainActor->VisibilityOn();
    m_wireframeActor->VisibilityOn();
    m_wireframeActor->GetProperty()->SetLineWidth( 2*line_w );
    break;
  }
  emit ActorUpdated();
}

void LayerSurface::SetActiveAnnotation( int n )
{
  if ( n < (int)m_annotations.size() )
  {
    if ( m_nActiveAnnotation < 0 && n >= 0 && this->GetProperty()->GetCurvatureMap() == LayerPropertySurface::CM_Threshold)
    {
      this->GetProperty()->SetCurvatureMap( LayerPropertySurface::CM_Binary );
    }
    m_nActiveAnnotation = n;
    SurfaceAnnotation* annot = GetActiveAnnotation();
    if (annot)
    {
      MRIS* mris = GetSourceSurface()->GetMRIS();
      if (mris->ct)
        CTABfree(&mris->ct);
      mris->ct = CTABdeepCopy(annot->GetColorTable());
      int* data = annot->GetAnnotationData();
      for (int i = 0; i < mris->nvertices; i++)
        mris->vertices[i].annotation = data[i];
    }
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
  double ratio = 1;
#if VTK_MAJOR_VERSION > 7
  ratio = MainWindow::GetMainWindow()->devicePixelRatio();
#endif
  m_vertexActor->SetVisibility( GetProperty()->GetShowVertices()? 1: 0 );
  m_vertexActor->GetProperty()->SetPointSize( GetProperty()->GetVertexPointSize()*ratio );
  m_vertexActor->GetProperty()->SetColor( GetProperty()->GetVertexColor() );
  for (int i = 0; i < 3; i++)
  {
    m_vertexActor2D[i]->SetVisibility( GetProperty()->GetShowVertices()? 1: 0 );
    m_vertexActor2D[i]->GetProperty()->SetPointSize( GetProperty()->GetVertexPointSize()*ratio );
    m_vertexActor2D[i]->GetProperty()->SetColor( GetProperty()->GetVertexColor() );
    if (GetProperty()->GetShowVertices())
        OnSlicePositionChanged(i);
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
  for (int i = 0; i < m_paths.size(); i++)
  {
    m_paths[i]->GetActor()->SetPosition(pos);
  }
  if (m_marks)
    m_marks->GetActor()->SetPosition(pos);

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

void LayerSurface::MoveLabelUp(SurfaceLabel *label)
{
  SurfaceLabel* activeLabel = GetActiveLabel();
  for (int i = 0; i < m_labels.size(); i++)
  {
    if (label == m_labels[i] && i > 0)
    {
      m_labels.removeAt(i);
      m_labels.insert(i-1, label);
      SetActiveLabel(activeLabel);

      UpdateOverlay(false);
      emit Modified();
      emit ActorChanged();
      return;
    }
  }
}

void LayerSurface::MoveLabelDown(SurfaceLabel *label)
{
  SurfaceLabel* activeLabel = GetActiveLabel();
  for (int i = 0; i < m_labels.size(); i++)
  {
    if (label == m_labels[i] && i < m_labels.size()-1)
    {
      m_labels.removeAt(i);
      m_labels.insert(i+1, label);
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
      label->deleteLater();
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

void LayerSurface::SetActiveSpline(SurfaceSpline *spline)
{
  for (int i = 0; i < m_splines.size(); i++)
  {
    if (spline == m_splines[i])
    {
      SetActiveSpline(i);
      return;
    }
  }
}


void LayerSurface::SetActiveSpline(int n)
{
  if (n >= 0 && n < m_splines.size() && n != m_nActiveSpline)
  {
    m_nActiveSpline = n;
    emit ActiveSplineChanged( n );
  }
}

void LayerSurface::MapLabels( unsigned char* data, int nVertexCount )
{
  for ( int i = m_labels.size()-1; i >= 0; i-- )
  {
    if (m_labels[i]->IsVisible())
      m_labels[i]->MapLabel( data, nVertexCount );
  }

  for ( int i = m_mappedLabels.size()-1; i >= 0; i-- )
  {
    if (m_mappedLabels[i]->IsVisible())
      m_mappedLabels[i]->MapLabelColorData( data, nVertexCount );
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

SurfaceSpline* LayerSurface::GetSpline(int n)
{
  if (n >= 0 && n < m_splines.size())
    return m_splines[n];
  else
    return NULL;
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
  return (GetFileName().toLower().contains("inflated") || GetActiveSurface() == FSSurface::SurfaceInflated);
}

bool LayerSurface::GetActiveLabelCentroidPosition(double *pos)
{
  SurfaceLabel* label = GetActiveLabel();
  int vno;
  double x, y, z;
  if (label && label->GetCentroid(&x, &y, &z, &vno))
  {
    return GetTargetAtVertex(vno, pos);
  }
  return false;
}

void LayerSurface::RemoveCurrentOverlay()
{
  if (m_nActiveOverlay >= 0)
  {
    m_overlays.removeAt(m_nActiveOverlay);
    SetActiveOverlay(m_overlays.size()-1);
  }
}

void LayerSurface::SetVisibleIn3D(bool bVisible)
{
  if (bVisible != m_bVisibleIn3D)
  {
    m_bVisibleIn3D = bVisible;
    emit VisibilityChanged(bVisible);
    emit ActorChanged();
  }
}

void LayerSurface::GetSmoothedVertexNormal(int nVertex, double *v_out)
{
  if (nVertex >= 0 && nVertex < this->GetNumberOfVertices())
    m_surfaceSource->GetSmoothedNormal(nVertex, v_out);
}

bool LayerSurface::LoadRGBFromFile(const QString &filename)
{
  QString fn = filename;
  fn.replace("~", QDir::homePath());
  if (!QFile::exists(fn))
  {
    fn = QFileInfo(QFileInfo(m_sFilename).dir(), filename).absoluteFilePath();
  }

  RGBMap map;
  map.name = QFileInfo(filename).completeBaseName();
  if (QFileInfo(fn).suffix() == "txt")
  {
    QFile file(fn);
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
      return false;
    while (!file.atEnd())
    {
      QString line = file.readLine();
      QStringList list = line.split(",", QString::SkipEmptyParts);
      if (list.size() < 3)
        list = line.split(" ", QString::SkipEmptyParts);
      if (list.size() == 3)
      {
        for (int i = 0; i < 3; i++)
          map.data << (int)list[i].toDouble();
      }
    }
    if (map.data.size() != GetNumberOfVertices()*3)
    {
      cout << "data size does not match" << endl;
      return false;
    }
  }
  else
  {
    MRI* mri = ::MRIread( filename.toLatin1().data() );
    if (!mri)
      return false;
    else if (mri->width != GetNumberOfVertices() || mri->height != 3)
    {
      cout << "data size does not match" << endl;
      MRIfree(&mri);
      return false;
    }

    switch ( mri->type )
    {
    case MRI_UCHAR:
      for (int i = 0; i < GetNumberOfVertices(); i++)
        for (int j = 0; j < 3; j++)
          map.data << MRIseq_vox( mri, i, j, 0, 0 );
      break;
    case MRI_INT:
      for (int i = 0; i < GetNumberOfVertices(); i++)
        for (int j = 0; j < 3; j++)
          map.data << MRIIseq_vox( mri, i, j, 0, 0 );
      break;

    case MRI_LONG:
      for (int i = 0; i < GetNumberOfVertices(); i++)
        for (int j = 0; j < 3; j++)
          map.data << MRILseq_vox( mri, i, j, 0, 0 );
      break;

    case MRI_FLOAT:
      for (int i = 0; i < GetNumberOfVertices(); i++)
        for (int j = 0; j < 3; j++)
          map.data << (int)MRIFseq_vox( mri, i, j, 0, 0 );
      break;

    case MRI_SHORT:
      for (int i = 0; i < GetNumberOfVertices(); i++)
        for (int j = 0; j < 3; j++)
          map.data << MRISseq_vox( mri, i, j, 0, 0 );
      break;
    default:
      MRIfree(&mri);
      return false;
    }
    MRIfree(&mri);
  }
  m_rgbMaps << map;
  SetActiveRGBMap(m_rgbMaps.size()-1);
  emit SurfaceRGBAdded();
  emit Modified();
  return true;
}

void LayerSurface::SetActiveRGBMap(int n)
{
  if ( n < m_rgbMaps.size() )
  {
    if ( m_nActiveRGBMap < 0 && n >= 0 )
    {
      this->GetProperty()->blockSignals(true);
      if (this->GetProperty()->GetCurvatureMap() == LayerPropertySurface::CM_Threshold)
        this->GetProperty()->SetCurvatureMap( LayerPropertySurface::CM_Binary );
      this->GetProperty()->blockSignals(false);
    }
    m_nActiveRGBMap = n;
    UpdateOverlay(false);
    emit ActiveOverlayChanged( m_nActiveOverlay );
    emit ActorUpdated();
    //      GetProperty()->SetShowOverlay(true);
    emit RGBMapChanged();
  }
}

QStringList LayerSurface::GetRGBMapNames()
{
  QStringList list;
  for (int i = 0; i < m_rgbMaps.size(); i++)
    list << m_rgbMaps[i].name;
  return list;
}

void LayerSurface::AddMappedLabel(LayerROI *label)
{
  if (!m_mappedLabels.contains(label))
  {
    m_mappedLabels << label;
    connect(label, SIGNAL(destroyed(QObject*)), SLOT(RemoveMappedLabel(QObject*)), Qt::UniqueConnection);
    connect(label, SIGNAL(VisibilityChanged(bool)), SLOT(UpdateOverlayLabels()), Qt::UniqueConnection);
    connect(label->GetProperty(), SIGNAL(ColorMapChanged()), SLOT(UpdateOverlayLabels()), Qt::UniqueConnection);
    connect(label->GetProperty(), SIGNAL(ThresholdChanged(double)), SLOT(UpdateOverlayLabels()), Qt::UniqueConnection);
    connect(this, SIGNAL(destroyed(QObject*)), label, SLOT(OnSurfaceDestroyed(QObject*)));
  }
}

void LayerSurface::RemoveMappedLabel(QObject *label_in)
{
  LayerROI* label = qobject_cast<LayerROI*>(label_in);
  if (label)
  {
    if (m_mappedLabels.contains(label))
    {
      m_mappedLabels.removeAll(label);
      disconnect(label, 0, this, 0);
    }
  }
  else
  {
    for (int i = 0; i < m_mappedLabels.size(); i++)
    {
      if (m_mappedLabels[i] == sender())
      {
        m_mappedLabels.removeAt(i);
        i--;
      }
    }
  }
  UpdateOverlay(true, true);
}

QVector<int> LayerSurface::FindPath(const QVector<int>& seeds)
{
  int* vert_vno = new int[seeds.size()];
  for (int i = 0; i < seeds.size(); i++)
    vert_vno[i] = seeds[i];

  int path[1000], path_length = 0;
  QVector<int> out_vno;
  if (m_surfaceSource->FindPath(vert_vno, seeds.size(), path, &path_length))
  {
    for (int i = 0; i < path_length; i++)
      out_vno << path[i];
  }
  delete[] vert_vno;
  return out_vno;
}

void LayerSurface::SetNeighborhoodSize(int nSize)
{
  MRIS* mris = m_surfaceSource->GetMRIS();
  if (mris->nsize == nSize)
    return;

  ::MRISsetNeighborhoodSizeAndDist(mris, nSize);
}

QVector<int> LayerSurface::GetVertexNeighbors(int vno)
{
  QVector<int> vno_list;
  MRIS* mris = m_surfaceSource->GetMRIS();
  VERTEX_TOPOLOGY* vt = &mris->vertices_topology[vno];
  for (int i = 0; i < vt->vtotal; i++)
    vno_list << vt->v[i];
  return vno_list;
}

void LayerSurface::ResetContralateralInfo()
{
  m_surfaceContralateral = NULL;
  m_surfaceSphere1 = NULL;
  m_surfaceSphere2 = NULL;
}

void LayerSurface::SetContralateralLayer(LayerSurface* layer, LayerSurface* sphere1, LayerSurface* sphere2)
{
  m_surfaceContralateral = layer;
  if (GetHemisphere() == sphere1->GetHemisphere())
  {
    m_surfaceSphere1 = sphere1;
    m_surfaceSphere2 = sphere2;
  }
  else
  {
    m_surfaceSphere1 = sphere2;
    m_surfaceSphere2 = sphere1;
  }
}

int LayerSurface::GetContralateralVertex(int vno)
{
  if (vno >= 0 && m_surfaceSphere1)
  {
    double ras[3];
    m_surfaceSphere1->GetSurfaceRASAtVertex(vno, ras);
    vno = m_surfaceSphere2->GetVertexAtSurfaceRAS(ras, NULL);
    return vno;
  }
  return vno;
}

bool LayerSurface::IsContralateralPossible()
{
  if (IsContralateralReady())
    return true;

  QString fn = GetFileName();
  QString fullpath = QFileInfo(fn).absolutePath();
  if (GetHemisphere() == 0)
    fn.replace("lh.", "rh.");
  else
    fn.replace("rh.", "lh.");

  return QFile::exists(fn) && QFile::exists(fullpath + "/lh.sphere.d1.left_right") &&
      QFile::exists(fullpath + "/rh.sphere.d1.left_right");
}

void LayerSurface::AddPathPoint(int vno)
{
  EditPathPoint(vno, false);
}

void LayerSurface::RemovePathPoint(int vno)
{
  EditPathPoint(vno, true);
}

void LayerSurface::ClearMarks()
{
  if (m_marks)
  {
    m_marks->Clear();
    emit ActorChanged();
  }
}

void LayerSurface::EditPathPoint(int vno, bool remove)
{
  if (!m_marks)
  {
    m_marks = new SurfacePath(this);
    m_marks->SetColor(Qt::cyan);
    //    SetActivePath(m_paths.size()-1);
    connect(m_marks, SIGNAL(Updated()), this, SIGNAL(ActorUpdated()));
    connect(m_marks, SIGNAL(CutLineMade()), this, SLOT(OnPathCut()), Qt::QueuedConnection);
    connect(m_marks, SIGNAL(PathMade()), this, SLOT(OnPathMade()), Qt::QueuedConnection);
    emit ActorChanged();
  }
  if (remove)
    m_marks->RemovePoint(vno);
  else
    m_marks->AddPoint(vno);
}

void LayerSurface::RemoveLastPathPoint()
{
  if (m_marks)
    m_marks->RemoveLastPoint();
}

void LayerSurface::SetActivePath(int n)
{
  if (m_nActivePath >=0)
    m_paths[m_nActivePath]->SetColor(Qt::red);

  if (n >= 0)
  {
    m_nActivePath = n;
    m_paths[n]->SetColor(Qt::yellow);
  }
  emit ActorUpdated();
}

SurfacePath* LayerSurface::GetActivePath()
{
  if (m_nActivePath >= 0)
    return m_paths[m_nActivePath];
  else
    return NULL;
}

SurfacePath* LayerSurface::GetMadePath(int nPath)
{
    if (nPath >= 0 && nPath < m_paths.size() && m_paths[nPath]->IsPathMade())
        return m_paths[nPath];
    else
        return NULL;
}

void LayerSurface::DeleteActivePath()
{
  if (m_nActivePath >= 0)
  {
    m_paths[m_nActivePath]->deleteLater();
    m_paths.removeAt(m_nActivePath);
    m_nActivePath = -1;
    emit ActorChanged();
  }
}

int LayerSurface::FindPathAt(int vno)
{
  QVector<int> verts = GetVertexNeighbors(vno);
  for (int i = 0; i < m_paths.size(); i++)
  {
    for (int j = 0; j < verts.size(); j++)
      if (m_paths[i]->Contains(verts[j]))
        return i;
  }
  return -1;
}

void LayerSurface::OnPathCut()
{
  SurfacePath* path = qobject_cast<SurfacePath*>(sender());
  if (path)
  {
    QVector<int> undoableVerts = m_surfaceSource->MakeCutLine(path->GetPathVerts());
    m_surfaceSource->RipFaces();
    m_surfaceSource->UpdatePolyData();
    path->SetUndoVerts(undoableVerts);
    emit ActorUpdated();
    PushMarksToPath();
  }
}

void LayerSurface::OnPathMade()
{
  PushMarksToPath();
}

void LayerSurface::PushMarksToPath()
{
  if (m_marks)
  {
    //    m_marks->SetColor(Qt::yellow);
    m_paths << m_marks;
    m_marks = NULL;
    SetActivePath(m_paths.size()-1);
  }
}

void LayerSurface::ClearAllCuts()
{
  //  for (int i = 0; i < m_paths.size(); i++)
  //  {
  //    m_paths[i]->deleteLater();
  //  }
  //  m_paths.clear();
  //  m_nActivePath = -1;
  m_surfaceSource->ClearCuts();
  m_surfaceSource->RipFaces();
  m_surfaceSource->UpdateHashTable();
  m_surfaceSource->UpdatePolyData();
  emit ActorChanged();
}

bool LayerSurface::HasUndoableCut()
{
  for (int i = 0; i < m_paths.size(); i++)
  {
    if (m_paths[i]->IsCutLineMade())
      return true;
  }
  return false;
}

void LayerSurface::UndoCut()
{
  for (int i = m_paths.size()-1; i >= 0; i--)
  {
    if (m_paths[i]->IsCutLineMade())
    {
      m_surfaceSource->ClearCuts(m_paths[i]->GetUndoVerts());
      m_paths[i]->deleteLater();
      m_paths.removeAt(i);
      if (m_nActivePath >= m_paths.size())
        m_nActivePath = -1;
      m_surfaceSource->RipFaces();
      m_surfaceSource->UpdateHashTable();
      m_surfaceSource->UpdatePolyData();
      emit ActorUpdated();
      return;
    }
  }
}

QVector<int> LayerSurface::FloodFillFromSeed(int seed_vno, const QVariantMap& options)
{
  MRIS* mris = m_surfaceSource->GetMRIS();
  char* filled;
  int num_filled_this_iter;
  int num_filled;
  int iter;
  int min_vno, max_vno, step_vno;
  int vno;
  //  int this_label = 0;
  int neighbor_index;
  int neighbor_vno;
  VERTEX* v;
  VERTEX_TOPOLOGY* vt;
  VERTEX* neighbor_v;
  //  float fvalue = 0;
  //  float seed_curv = 0;
  //  float seed_fvalue = 0;
  //  int new_index;
  //  int num_labels_found, found_label_index;
  //  int skip;
  int count;

  QVector<int> filled_verts;

  if (seed_vno < 0 || seed_vno >= mris->nvertices)
    return filled_verts;

  /* init filled array. */
  filled = (char*) calloc (mris->nvertices, sizeof(char));
  memset(filled, 0, sizeof(char)*mris->nvertices);

  /* start with the seed filled.*/
  filled[seed_vno] = TRUE;
  filled_verts << seed_vno;

  /* find seed values for some conditions. */
  //  if (params->dont_cross_label)
  //    this_label = labl_selected_label;
  //  if (params->dont_cross_cmid)
  //    seed_curv = mris->vertices[seed_vno].curv;
  //  if (params->dont_cross_fthresh)
  //    sclv_get_value (&mris->vertices[seed_vno],
  //                    sclv_current_field, &seed_fvalue);

  /* while we're still filling stuff in a pass... */
  num_filled_this_iter = 1;
  num_filled = 0;
  iter = 0;
  bool bDoNotCrossPaths = options["DoNotCrossPaths"].toBool();
  bool bDoNotCrossLabels = options["DoNotCrossLabels"].toBool();
  bool bDoNotFillUnlabeled = options["DoNotFillUnlabeled"].toBool();
  bool bDoNotCrossThreshold = options["DoNotCrossThreshold"].toBool();
  bool bDoNotCrossCurv = (HasCurvature() && options["FillToCurvature"].toBool());
  bool bNewLabel = options["CreateLabel"].toBool();
  bool bAsAnnotation = options["AsAnnotation"].toBool();
  double seed_curv = 0;
  double seed_value = 0;
  SurfaceOverlay* overlay = GetActiveOverlay();
  if (!overlay)
    bDoNotCrossThreshold = false;
  double fthresh = 0;
  if (bDoNotCrossThreshold && overlay)
  {
    seed_value = overlay->GetDataAtVertex(seed_vno);
    fthresh = overlay->GetProperty()->GetMinPoint();
  }
  if (bDoNotCrossCurv)
    seed_curv = mris->vertices[seed_vno].curv;

  while (num_filled_this_iter > 0)
  {
    num_filled_this_iter = 0;

    /* switch between iterating forward and backwards. */
    if ((iter%2)==0)
    {
      min_vno = 0;
      max_vno = mris->nvertices-1;
      step_vno = 1;
    }
    else
    {
      min_vno = mris->nvertices-1;
      max_vno = 0;
      step_vno = -1;
    }

    /* for each vertex, if it's filled, check its neighbors. for the
       rules that are up-to-and-including, make the check on this
       vertex. for the rules that are up-to-and-not-including, check
       on the neighbor. */
    for (vno = min_vno; vno != max_vno; vno += step_vno)
    {
      if (filled[vno])
      {

        /* check the neighbors... */
        v = &mris->vertices[vno];
        vt = &mris->vertices_topology[vno];

        /* if this vert is ripped, move on. */
        if (v->ripflag)
        {
          continue;
        }

        /* if we're not crossing paths, check if this is a
           path. if so, move on. */
        if (bDoNotCrossPaths && IsVertexOnPath(vno))
        {
          continue;
        }

        if (bDoNotCrossCurv && ((seed_curv <= 0 && v->curv > 0) || (seed_curv >= 0 && v->curv < 0)))
          continue;

        /* if we're not crossing the cmid, see if the cmid at this
           vertex is on the other side of the cmid as the seed
           point. if so, move on. */
        //        if (params->dont_cross_cmid &&
        //            ((seed_curv <= cmid && v->curv > cmid) ||
        //             (seed_curv >= cmid && v->curv < cmid)))
        //        {
        //          continue;
        //        }

        for (neighbor_index = 0;
             neighbor_index < vt->vnum;
             neighbor_index++)
        {
          neighbor_vno = vt->v[neighbor_index];
          neighbor_v = &mris->vertices[neighbor_vno] ;

          /* if the neighbor is filled, move on. */
          if (filled[neighbor_vno])
            continue;

          if (neighbor_v->ripflag)
            continue;

          /* if we're not crossing labels, check if the label at
             this vertex is the same as the one at the seed. if not,
             move on. */
          if (bDoNotCrossLabels || bDoNotFillUnlabeled)
          {
            bool bFound = false, bSkip = false;
            if (bAsAnnotation)
            {
              SurfaceAnnotation* annot = GetActiveAnnotation();
              int nIndex = annot->GetIndexAtVertex(neighbor_vno);
              if (nIndex >= 0)
              {
                bFound = true;
                if (nIndex != annot->GetIndexAtVertex(seed_vno))
                {
                  bSkip = true;
                }
              }
            }
            else
            {
              for (int i = 0; i < m_labels.size(); i++)
              {
                if (m_labels[i]->HasVertex(neighbor_vno))
                {
                  bFound = true;
                  if (i != m_nActiveLabel || bNewLabel)
                  {
                    bSkip = true;
                    break;
                  }
                }
              }
            }

            if (bSkip && bDoNotCrossLabels)
              continue;

            if (!bFound && bDoNotFillUnlabeled)
              continue;
          }


          /* if we're not crossing the fthresh, make sure this
             point is above it, or, if our initial functional
             value was negative, make sure it's not above
             -fthresh. if not, move on. */
          if (bDoNotCrossThreshold)
          {
            double fvalue = overlay->GetDataAtVertex(neighbor_vno);
            if ((fthresh != 0 &&
                 seed_value > 0 &&
                 fvalue < fthresh) ||
                (fthresh != 0 &&
                 seed_value < 0 &&
                 fvalue > -fthresh) ||
                (fthresh == 0 && (fvalue * seed_value < 0)))
            {
              continue;
            }
          }

          /* mark this vertex as filled. */
          filled[neighbor_vno] = TRUE;
          filled_verts << neighbor_vno;
          num_filled_this_iter++;
          num_filled++;
        }
      }
    }

    iter++;
  }

  /* mark all filled vertices. */
  if (!options["DoNotMarkSurface"].toBool())
  {
    for (vno = 0; vno < mris->nvertices; vno++ )
    {
      mris->vertices[vno].ripflag = (!filled[vno]);
    }
  }

  free (filled);

  if (filled_verts.size() == 1)
    filled_verts.clear();

  return filled_verts;
}

bool LayerSurface::FillUncutArea(int vno)
{
  if (vno < 0)
    return false;

  FloodFillFromSeed(vno);
  m_surfaceSource->RipFaces();
  m_surfaceSource->UpdateHashTable();
  m_surfaceSource->UpdatePolyData();

  if (m_nActivePath >=0 && !m_paths[m_nActivePath]->IsCutLineMade())
  {
    m_paths[m_nActivePath]->deleteLater();
    m_paths.removeAt(m_nActivePath);
    m_nActivePath = -1;
    emit ActorChanged();
  }
  else
    emit ActorUpdated();
  return true;
}

bool LayerSurface::LoadPatch(const QString &filename)
{
  if (m_surfaceSource->LoadPatch(filename))
  {
    m_sPatchFilename = filename;
    emit ActorUpdated();
    return true;
  }
  else
    return false;
}

bool LayerSurface::WritePatch(const QString &filename)
{
  return (::MRISwritePatch(m_surfaceSource->GetMRIS(), filename.toLatin1().data()) == 0);
}

bool LayerSurface::IsVertexRipped(int vno)
{
  for (int i = 0; i < m_paths.size(); i++)
  {
    if (m_paths[i]->IsCutLineMade() && m_paths[i]->Contains(vno))
      return true;
  }
  return false;
}

bool LayerSurface::IsVertexOnPath(int vno)
{
  foreach (SurfacePath* path, m_paths)
  {
    if (path->IsPathMade() && path->Contains(vno))
      return true;
  }
  return false;
}

bool LayerSurface::FillPath(int nvo, const QVariantMap &options)
{
  bool bAsAnnotation = options["AsAnnotation"].toBool();
  QVector<int> verts = FloodFillFromSeed(nvo, options);
  if (verts.size() == 0)
  {
    cout << "Did not fill/remove any vertices" << endl;
    return false;
  }

  if (bAsAnnotation)
  {
    SurfaceAnnotation* annot = GetActiveAnnotation();
    annot->EditLabel(verts, options["FillAnnotationIndex"].toInt(), options);
    emit Modified();
  }
  else
  {
    if (options["CreateLabel"].toBool())
      CreateNewLabel();

    SurfaceLabel* label = GetActiveLabel();
    if (label)
    {
      //    label->SaveForUndo();
      label->EditVertices(verts, !options["RemoveFromLabel"].toBool());
      emit Modified();
    }
  }
  return true;
}

int LayerSurface::FillPath(const QVector<int> &verts, const QVariantMap &options)
{
  QVariantMap opt = options;
  foreach (int nvo, verts)
  {
    if (FillPath(nvo, opt) && opt["CreateLabel"].toBool())
    {
      opt["CreateLabel"] = false;
      opt["AddToLabel"] = true;
    }
  }
  return opt.value("FillAnnotationIndex").toInt();
}

int LayerSurface::GetLastMark()
{
  int vno = -1;
  if (m_marks)
  {
    QVector<int> verts = m_marks->GetPathVerts();
    if (!verts.isEmpty())
      vno = verts.last();
  }
  return vno;
}

QVector<int> LayerSurface::GetAllMarks()
{
  QVector<int> verts;
  if (m_marks)
    verts = m_marks->GetPathVerts();
  return verts;
}

bool LayerSurface::LoadParameterization(const QString &filename)
{
  MRIS* mris = m_surfaceSource->GetMRIS();
  MRI* mri = ::MRISreadParameterizationToSurface(mris, filename.toLatin1().data() );
  if ( mri )
  {
    int nPerFrame = mri->width*mri->height*mri->depth;
    int nframes = nPerFrame*mri->nframes / mris->nvertices;
    float* data = new float[nPerFrame*mri->nframes];
    for (int nx = 0; nx < mri->width; nx++)
    {
      for (int ny = 0; ny < mri->height; ny++)
      {
        for (int nz = 0; nz < mri->depth; nz++)
        {
          for (int nk = 0; nk < mri->nframes; nk++)
          {
            data[nk*nPerFrame + nz*mri->height*mri->width + ny*mri->width + nx]
                = ::MRIgetVoxVal(mri, nx, ny, nz, nk);
          }
        }
      }
    }
    // create overlay
    SurfaceOverlay* overlay = new SurfaceOverlay( this );
    overlay->InitializeData(data, mris->nvertices, nframes);
    overlay->SetName( QFileInfo(filename).fileName() );
    overlay->SetFileName( filename );

    m_overlays.push_back( overlay );
    SetActiveOverlay( m_overlays.size() - 1 );

    emit Modified();
    emit SurfaceOverlayAdded( overlay );
    connect(overlay, SIGNAL(DataUpdated()), this, SIGNAL(SurfaceOverlyDataUpdated()), Qt::UniqueConnection);
    MRIfree(&mri);
    return true;
  }
  else
    return false;
}

bool LayerSurface::LoadCoordsFromParameterization(const QString &filename)
{
  MRIS* mris = m_surfaceSource->GetMRIS();
  MRI_SP* mrisp = ::MRISPread(filename.toLatin1().data());
  if (mrisp)
  {
    ::MRIScoordsFromParameterization(mrisp, mris, IsInflated()?WHITE_VERTICES:CURRENT_VERTICES);
    m_surfaceSource->UpdateHashTable();
    m_surfaceSource->UpdateCoords();
    emit ActorUpdated();
    return true;
  }
  else
  {
    cerr << "Failed to load " << qUtf8Printable(filename) << endl;
    return false;
  }
}

vtkActor* LayerSurface::GetMainActor()
{
  return m_mainActor;
}

void LayerSurface::SetHighlightedLabelOnAnnotation(int n)
{
  if (m_nActiveAnnotation >= 0)
    GetActiveAnnotation()->SetHighlightedLabel(n);
}

SurfaceAnnotation* LayerSurface::CreateNewAnnotation(const QString &ct_file, const QString& name)
{
  // create annotation
  SurfaceAnnotation* annot = new SurfaceAnnotation( this );
  if (annot->InitializeNewAnnotation(ct_file))
  {
    annot->SetName(name.isEmpty()?"Unnamed":name);
    m_annotations.push_back( annot );
    connect(annot, SIGNAL(Modified()), SLOT(UpdateColorMap()));
    SetActiveAnnotation( m_annotations.size() - 1 );

    emit Modified();
    emit SurfaceAnnotationAdded( annot );
    emit ActorUpdated();
  }
  else
  {
    delete annot;
    annot = NULL;
  }
  return annot;
}
