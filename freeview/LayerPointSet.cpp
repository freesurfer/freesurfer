/**
 * @brief Layer data object for MRI surface.
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

#include "LayerPointSet.h"
#include "LayerMRI.h"
#include "FSPointSet.h"
#include "LayerPropertyPointSet.h"
#include "FSVolume.h"
#include "vtkRGBAColorTransferFunction.h"
#include <vtkRenderer.h>
#include <vtkActor.h>
#include <vtkMath.h>
#include <vtkSphereSource.h>
#include <vtkSplineFilter.h>
#include <vtkPoints.h>
#include <vtkPolyDataMapper.h>
#include <vtkAppendPolyData.h>
#include <vtkCellArray.h>
#include <vtkFloatArray.h>
#include <vtkTubeFilter.h>
#include <vtkPolyData.h>
#include <vtkProperty.h>
#include <vtkPlane.h>
#include <vtkCutter.h>
#include <vtkStripper.h>
#include <vtkTriangleFilter.h>
#include <vtkSmartPointer.h>
#include <vtkPointData.h>
#include <QDebug>
#include <QFile>
#include <QJsonDocument>
#include "MyUtils.h"
#include "vtkContourTriangulator.h"

#define NUM_OF_SIDES  10  // must be even number!
#define POINTSET_VERSION  1

LayerPointSet::LayerPointSet( LayerMRI* ref, int nType, QObject* parent ) : LayerEditable( parent )
{
  m_strTypeNames.push_back( "PointSet" );
  m_sPrimaryType = "PointSet";
  m_actorBalls = vtkActor::New();
  m_actorSpline = vtkActor::New();
  for ( int i = 0; i < 3; i++ )
  {
    m_actorSlice[i] = vtkActor::New();
    m_actorSlice[i]->GetProperty()->SetInterpolationToFlat();
    m_actorSlice[i]->GetProperty()->SetAmbient( 1 );
    m_actorSlice[i]->GetProperty()->SetDiffuse( 0 );
    m_actorSplineSlice[i] = vtkActor::New();
    m_actorSplineSlice[i]->GetProperty()->SetInterpolationToFlat();
    m_actorSplineSlice[i]->GetProperty()->SetAmbient( 1 );
    m_actorSplineSlice[i]->GetProperty()->SetDiffuse( 0 );
    double pos[3] = {0,0,0};
    pos[i] = 1e-4;
    if (i == 2)
      pos[i] = -pos[i];
    m_actorSplineSlice[i]->SetPosition(pos);
  }
  m_layerRef = ref;
  m_pointSetSource = new FSPointSet();

  m_splinedPoints = vtkSmartPointer<vtkPoints>::New();

  mProperty = new LayerPropertyPointSet( this );
  GetProperty()->SetType( nType );

  m_mapper = vtkSmartPointer<vtkPolyDataMapper>::New();

  if ( nType == LayerPropertyPointSet::ControlPoint || nType == LayerPropertyPointSet::Enhanced)
  {
    GetProperty()->SetShowSpline( false );
    GetProperty()->SetRadius ( 0.5 );
    GetProperty()->SetSnapToVoxelCenter( true );
    GetProperty()->SetColor( 0, 1, 0 );
  }

  LayerPropertyPointSet* p = GetProperty();
  connect(p, SIGNAL(OpacityChanged(double)), this, SLOT(UpdateOpacity()));
  connect(p, SIGNAL(ColorMapChanged()), this, SLOT(UpdateColorMap()));
  connect(p, SIGNAL(RadiusChanged(double)), this, SLOT(RebuildActors()));
  connect(p, SIGNAL(SplineRadiusChanged(double)), this, SLOT(RebuildActors()));
  connect(p, SIGNAL(ScalarChanged()), this, SLOT(RebuildActors()));
  connect(p, SIGNAL(SnapToVoxelCenterChanged(bool)), this, SLOT(UpdateSnapToVoxelCenter()));
  connect(p, SIGNAL(SplineVisibilityChanged(bool)), this, SLOT(UpdateSplineVisibility()));
  connect(p, SIGNAL(ClosedSplineChanged(bool)), this, SLOT(RebuildActors()));
}

LayerPointSet::~LayerPointSet()
{
  m_actorBalls->Delete();
  m_actorSpline->Delete();
  for ( int i = 0; i < 3; i++ )
  {
    m_actorSlice[i]->Delete();
    m_actorSplineSlice[i]->Delete();
  }
}

bool LayerPointSet::LoadFromFile( const QString& filename )
{
  if (!LoadFromJsonFile(filename))
  {
    if ( GetProperty()->GetType() == LayerPropertyPointSet::ControlPoint )
    {
      if ( !m_pointSetSource->ReadAsControlPoints( filename ) )
      {
        return false;
      }
      else
      {
        cout << "Warning: Coordinate of control points has been converted to realRAS in "
             << qPrintable(filename) << " and will be saved in that way."<< endl << endl;
      }
    }
    else
    {
      if ( !m_pointSetSource->ReadAsLabel( filename ) )
      {
        return false;
      }
    }

    GetProperty()->SetStatRange(m_pointSetSource->GetMinStat(), m_pointSetSource->GetMaxStat());
    m_points.clear();
    m_pointSetSource->LabelToPointSet( m_points, m_layerRef->GetSourceVolume() );
  }
  SetFileName( filename );
  RebuildActors();

  return true;
}

bool LayerPointSet::LoadFromJsonFile(const QString &filename)
{
  QFile file( filename );
  if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
  {
    cerr << qPrintable(file.errorString()) << "\n";;
    return false;
  }
  m_mapEnhancedData = QJsonDocument::fromJson(file.readAll()).toVariant().toMap();
  if (m_mapEnhancedData.isEmpty())
    return false;

  if (m_mapEnhancedData.value("data_type").toString() != "fs_pointset")
  {
    cerr << "Not a freesurfer pointset file\n";
    return false;
  }

  QVariantList list = m_mapEnhancedData.value("points").toList();
  QString coord_strg = m_mapEnhancedData.value("vox2ras").toString();
  FSVolume* ref_vol = m_layerRef->GetSourceVolume();
  m_points.clear();
  foreach (QVariant v, list)
  {
    QVariantMap map = v.toMap();
    QVariantMap coords = map["coordinates"].toMap();
    ControlPoint wp;
    wp.pt[0] = coords["x"].toDouble();
    wp.pt[1] = coords["y"].toDouble();
    wp.pt[2] = coords["z"].toDouble();
    wp.value = map["legacy_stat"].toDouble();
    wp.info = map;
    if ( coord_strg == "tkreg" )
    {
      ref_vol->TkRegToNativeRAS( wp.pt, wp.pt );
    }
    else if (coord_strg == "voxel")
    {
      MRIvoxelToWorld(ref_vol->GetMRI(), wp.pt[0], wp.pt[1], wp.pt[2], wp.pt, wp.pt+1, wp.pt+2);
    }
    ref_vol->NativeRASToRAS( wp.pt, wp.pt );
    ref_vol->RASToTarget( wp.pt, wp.pt );
    m_points << wp;
  }
  GetProperty()->SetType(LayerPropertyPointSet::Enhanced);

  if (m_mapEnhancedData.contains("color"))
  {
    QVariantList list = m_mapEnhancedData["color"].toList();
    if (list.size() == 3)
    {
      GetProperty()->SetColor(QColor(list[0].toInt(), list[1].toInt(), list[2].toInt()));
    }
  }

  return true;
}

bool LayerPointSet::LoadFromString(const QString &content)
{
  if (!m_pointSetSource->ReadFromStringAsControlPoints(content))
    return false;

  m_points.clear();
  m_pointSetSource->LabelToPointSet( m_points, m_layerRef->GetSourceVolume() );
  RebuildActors();

  return true;
}

bool LayerPointSet::Save()
{
  if ( m_sFilename.size() == 0 || m_layerRef == NULL )
  {
    return false;
  }

  m_pointSetSource->UpdateLabel( m_points, m_layerRef->GetSourceVolume() );

  bool bSaved = false;
  if ( GetProperty()->GetType() == LayerPropertyPointSet::ControlPoint )
    bSaved = m_pointSetSource->WriteAsControlPoints( m_sFilename );
  else if (GetProperty()->GetType() == LayerPropertyPointSet::WayPoint )
    bSaved = m_pointSetSource->WriteAsLabel( m_sFilename );
  else
    bSaved = SaveAsJson(m_sFilename);

  if ( !bSaved )
  {
    m_bModified = true;
  }

  return bSaved;
}

bool LayerPointSet::SaveAsJson(const QString& filename)
{
  QVariantList list;
  FSVolume* ref_vol = m_layerRef->GetSourceVolume();
  double pos[3];
  foreach (ControlPoint p, m_points)
  {
    QVariantMap map = p.info;
    // convert to tkreg coords
    ref_vol->TargetToRAS( p.pt, pos );
    ref_vol->RASToNativeRAS( pos, pos );
//    ref_vol->NativeRASToTkReg(pos, pos);
    QVariantMap coords;
    coords["x"] = pos[0];
    coords["y"] = pos[1];
    coords["z"] = pos[2];
    map["coordinates"] = coords;
    if (!map.contains("legacy_stat"))
      map["legacy_stat"] = p.value;
    list << map;
  }

  m_mapEnhancedData["points"] = list;
  m_mapEnhancedData["data_type"] = "fs_pointset";
  m_mapEnhancedData["vox2ras"] = "scanner_ras";
  double* rgb = GetProperty()->GetColor();
  list.clear();
  QColor c = QColor::fromRgbF(rgb[0], rgb[1], rgb[2]);
  list << c.red() << c.green() << c.blue();
  m_mapEnhancedData["color"] = list;

  QFile file( filename );
  if (!file.open(QIODevice::WriteOnly | QIODevice::Text))
  {
    QString strg = file.errorString();
    if (strg.isEmpty())
      cerr << "Can not open file for writing\n";
    else
      cerr << qPrintable(strg) << "\n";
    return false;
  }

  m_mapEnhancedData["version"] = POINTSET_VERSION;
  file.write(QJsonDocument::fromVariant(m_mapEnhancedData).toJson());
  return true;
}

void LayerPointSet::UpdateLabelData()
{
  m_pointSetSource->UpdateLabel( m_points, m_layerRef->GetSourceVolume() );
}

bool LayerPointSet::HasUndo()
{
  return m_bufferUndo.size() > 0;
}

bool LayerPointSet::HasRedo()
{
  return m_bufferRedo.size() > 0;
}

void LayerPointSet::Undo()
{
  if ( m_bufferUndo.size() > 0 )
  {
    m_bufferRedo.push_back( m_points );
    m_points = m_bufferUndo[m_bufferUndo.size()-1];
    m_bufferUndo.pop_back();
    RebuildActors();
  }
}

void LayerPointSet::Redo()
{
  if ( m_bufferRedo.size() > 0 )
  {
    m_bufferUndo.push_back( m_points );
    m_points = m_bufferRedo[m_bufferRedo.size()-1];
    m_bufferRedo.pop_back();
    RebuildActors();
  }
}

void LayerPointSet::SaveForUndo()
{
  m_bufferUndo.push_back( m_points );
  m_bufferRedo.clear();
}

void LayerPointSet::Append2DProps( vtkRenderer* renderer, int nPlane )
{
  if ( nPlane < 3 && nPlane >= 0 )
  {
    renderer->AddViewProp( m_actorSplineSlice[nPlane] );
    renderer->AddViewProp( m_actorSlice[nPlane] );
  }
}

void LayerPointSet::Append3DProps( vtkRenderer* renderer, bool* bSliceVisibility )
{
  Q_UNUSED(renderer);
  Q_UNUSED(bSliceVisibility);
  renderer->AddViewProp( m_actorSpline );
  renderer->AddViewProp( m_actorBalls );
}

void LayerPointSet::SetVisible( bool bVisible )
{
  if ( GetProperty()->GetShowSpline() )
  {
    m_actorSpline->SetVisibility( bVisible ? 1 : 0 );
    for (int i = 0; i < 3; i++)
      m_actorSplineSlice[i]->SetVisibility( bVisible ? 1 : 0 );
  }
  else
  {
    m_actorSpline->SetVisibility( 0 );
    for (int i = 0; i < 3; i++)
      m_actorSplineSlice[i]->SetVisibility( 0 );
  }
  m_actorBalls->SetVisibility( bVisible ? 1 : 0 );
  for ( int i = 0; i < 3; i++ )
  {
    m_actorSlice[i]->SetVisibility( bVisible ? 1 : 0 );
  }

  LayerEditable::SetVisible(bVisible);
}

bool LayerPointSet::IsVisible()
{
  return m_actorSlice[0]->GetVisibility() > 0;
}

void LayerPointSet::UpdateSplineVisibility()
{
  SetVisible(IsVisible());
  RebuildActors();
}

bool LayerPointSet::HasProp( vtkProp* prop )
{
  return ( prop == m_actorSpline || prop == m_actorBalls );
}

void LayerPointSet::OnSlicePositionChanged( int nPlane )
{
  Q_UNUSED(nPlane);
  RebuildActors( false );   // no need to rebuild 3D points
}

void LayerPointSet::RebuildActors( bool bRebuild3D )
{
  if ( !m_layerRef )
  {
    return;
  }

  blockSignals( true );

  // 3D
  MRI* mri = m_layerRef->GetSourceVolume()->GetMRITarget();
  double voxel_size[3] = { mri->xsize, mri->ysize, mri->zsize };
  // double* origin = m_layerRef->GetWorldOrigin();
  double scale = qMin( voxel_size[0], qMin( voxel_size[1], voxel_size[2] ) );
  double radius = GetProperty()->GetRadius();
  bool bClosed = GetProperty()->GetClosedSpline();

  vtkSmartPointer<vtkAppendPolyData> append = vtkSmartPointer<vtkAppendPolyData>::New();
  vtkPoints* pts = vtkPoints::New();
  vtkCellArray* lines = vtkCellArray::New();
  lines->InsertNextCell( m_points.size() + (bClosed?1:0) );
  for ( int i = 0; i < m_points.size(); i++ )
  {
    if (radius > 0)
    {
      vtkSphereSource* sphere = vtkSphereSource::New();
      sphere->SetCenter( m_points[i].pt );
      sphere->SetRadius( radius * scale );
      sphere->SetThetaResolution( 10 );
      sphere->SetPhiResolution( 20 );
      append->AddInputConnection( sphere->GetOutputPort() );
      sphere->Delete();
    }
    pts->InsertNextPoint( m_points[i].pt);
    lines->InsertCellPoint( i );
  }
  if (bClosed && m_points.size() > 1)
  {
//    int nLast = m_points.size()-1;
//    // hack to avoid a VTK artifact
//    double dDist = sqrt(vtkMath::Distance2BetweenPoints( m_points[0].pt, m_points[nLast].pt ));
//    if (dDist == 0)
//      dDist = 1;
//    pts->InsertNextPoint( m_points[0].pt[0] - (m_points[0].pt[0] - m_points[nLast].pt[0])/dDist*0.1*scale,
//                          m_points[0].pt[1] - (m_points[0].pt[1] - m_points[nLast].pt[1])/dDist*0.1*scale,
//                          m_points[0].pt[2] - (m_points[0].pt[2] - m_points[nLast].pt[2])/dDist*0.1*scale );
//    lines->InsertCellPoint(m_points.size());
    lines->InsertCellPoint(0);
  }
  vtkPolyDataMapper* mapper = vtkPolyDataMapper::New();
  if ( m_points.size() > 0 && radius > 0 )
  {
    mapper->SetInputConnection( append->GetOutputPort() );
  }
  else
  {
#if VTK_MAJOR_VERSION > 5
    mapper->SetInputData( vtkSmartPointer<vtkPolyData>::New() );
#else
    mapper->SetInput( vtkSmartPointer<vtkPolyData>::New() );
#endif
  }
  m_actorBalls->SetMapper( mapper );
  mapper->Delete();

  vtkPolyData* polydata_tube = NULL;
  if ( GetProperty()->GetShowSpline() )
  {
    vtkPolyData* polydata = vtkPolyData::New();
    polydata->SetPoints( pts );
    polydata->SetLines( lines );
    vtkSplineFilter* spline = vtkSplineFilter::New();
//    spline->SetSubdivideToSpecified();
//    spline->SetNumberOfSubdivisions(pts->GetNumberOfPoints()*2);
    spline->SetSubdivideToLength();
    spline->SetLength(scale*2);
    if ( GetProperty()->GetScalarType() == LayerPropertyPointSet::ScalarSet )
    {
      spline->SetSubdivideToSpecified();
      spline->SetNumberOfSubdivisions( GetProperty()->GetActiveScalarSet().nNum - 1 );
    }
    if ( m_points.size() > 1 )
    {
    //  polydata->Update();
      UpdateScalars(polydata);
#if VTK_MAJOR_VERSION > 5
      spline->SetInputData( polydata );
#else
      spline->SetInput( polydata );
#endif
      vtkTubeFilter* tube = vtkTubeFilter::New();
      tube->SetNumberOfSides( NUM_OF_SIDES );
      tube->SetInputConnection( spline->GetOutputPort() );
      tube->SetRadius( GetProperty()->GetSplineRadius() * scale );
//      if (!bClosed)
        tube->CappingOn();
      m_mapper->SetInputConnection( tube->GetOutputPort() );
      tube->Update();
      m_splinedPoints = spline->GetOutput()->GetPoints();
      polydata_tube = tube->GetOutput();
      m_actorSpline->SetMapper( m_mapper );
      tube->Delete();
    }
    polydata->Delete();
    spline->Delete();
  }

  pts->Delete();
  lines->Delete();

  // 2D
  for ( int i = 0; i < 3; i++ )
  {
    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    vtkSmartPointer<vtkAppendPolyData> append = vtkSmartPointer<vtkAppendPolyData>::New();
    int n = 0;
    for ( int j = 0; j < m_points.size(); j++ )
    {
      if ( radius > 0 && fabs( m_dSlicePosition[i] - m_points[j].pt[i] ) < ( voxel_size[i] / 2 ) )
      {
        vtkSphereSource* sphere = vtkSphereSource::New();
        double point[3] = { m_points[j].pt[0], m_points[j].pt[1], m_points[j].pt[2] };
        point[i] = m_dSlicePosition[i];
        sphere->SetCenter( point );
        sphere->SetRadius( radius * scale );
        sphere->SetThetaResolution( 12 );
        sphere->SetPhiResolution( 24 );

        vtkSmartPointer<vtkPlane> plane = vtkSmartPointer<vtkPlane>::New();
        plane->SetOrigin( point );
        plane->SetNormal( i==0?1:0, i==1?1:0, i==2?1:0 );

        vtkSmartPointer<vtkCutter> cutter =
            vtkSmartPointer<vtkCutter>::New();
        cutter->SetInputConnection( sphere->GetOutputPort() );
        cutter->SetCutFunction( plane );

        vtkSmartPointer<vtkStripper> stripper = vtkSmartPointer<vtkStripper>::New();
        stripper->SetInputConnection( cutter->GetOutputPort() );
        stripper->Update();

        vtkSmartPointer<vtkPolyData> cutpoly = vtkSmartPointer<vtkPolyData>::New();
        cutpoly->SetPoints( stripper->GetOutput()->GetPoints() );
        cutpoly->SetPolys( stripper->GetOutput()->GetLines() );

        vtkSmartPointer<vtkTriangleFilter> triangleFilter = vtkSmartPointer<vtkTriangleFilter>::New();
#if VTK_MAJOR_VERSION > 5
        triangleFilter->SetInputData( cutpoly );
#else
        triangleFilter->SetInput( cutpoly );
#endif

        // append->AddInput( triangleFilter->GetOutput() );
        append->AddInputConnection( sphere->GetOutputPort() );
        sphere->Delete();
        n++;
      }
    }
    if ( n > 0 )
    {
      mapper->SetInputConnection( append->GetOutputPort() );
    }
    else
    {
#if VTK_MAJOR_VERSION > 5
      mapper->SetInputData( vtkSmartPointer<vtkPolyData>::New() );
#else
      mapper->SetInput( vtkSmartPointer<vtkPolyData>::New() );
#endif
    }
    m_actorSlice[i]->SetMapper( mapper );

    // 2D tube
    if (polydata_tube)
    {
      vtkSmartPointer<vtkPlane> plane = vtkSmartPointer<vtkPlane>::New();
      plane->SetOrigin( m_dSlicePosition );
      plane->SetNormal( i==0?1:0, i==1?1:0, i==2?1:0 );

      vtkSmartPointer<vtkCutter> cutter =
          vtkSmartPointer<vtkCutter>::New();
#if VTK_MAJOR_VERSION > 5
      cutter->SetInputData(polydata_tube);
#else
      cutter->SetInput(polydata_tube);
#endif
      cutter->SetCutFunction( plane );

      vtkSmartPointer<vtkStripper> stripper = vtkSmartPointer<vtkStripper>::New();
      stripper->SetInputConnection( cutter->GetOutputPort() );
      vtkSmartPointer<vtkContourTriangulator> ct = vtkSmartPointer<vtkContourTriangulator>::New();
      ct->SetInputConnection(stripper->GetOutputPort());
      mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
      mapper->SetInputConnection(ct->GetOutputPort());

      m_actorSplineSlice[i]->SetMapper(mapper);
    }
  }

  UpdateColorMap();

  if ( !bRebuild3D )
  {
    blockSignals(false);
    return;
  }

  UpdateOpacity();

  blockSignals(false);
  emit ActorUpdated();
}

void LayerPointSet::UpdateScalars(vtkPolyData* polydata)
{
  if ( true )
  {
    LayerMRI* layer = GetProperty()->GetScalarLayer();
    vtkPoints* pts = polydata->GetPoints();
    int nPts = pts->GetNumberOfPoints();
    vtkFloatArray* scalars = vtkFloatArray::New();
    scalars->SetNumberOfValues( nPts );
    //    double pt[3] = { 0, 0, 0 };
    double val = 0;
    int nNum = m_points.size();
    for ( int i = 0; i < m_points.size(); i++ )
    {
      if ( true ) // (i%NUM_OF_SIDES) == 0 )
      {
        //        double* p1 = pts->GetPoint( i );
        //        double* p2 = pts->GetPoint( i + NUM_OF_SIDES/2 );
        //        pt[0] = ( p1[0] + p2[0] ) / 2;
        //        pt[1] = ( p1[1] + p2[1] ) / 2;
        //        pt[2] = ( p1[2] + p2[2] ) / 2;
        if ( layer && GetProperty()->GetScalarType() == LayerPropertyPointSet::ScalarLayer )
        {
          val = layer->GetVoxelValue( pts->GetPoint(i) );
        }
        else if ( GetProperty()->GetScalarType() == LayerPropertyPointSet::ScalarSet )
        {
          val = GetProperty()->GetActiveScalarSet().dValue[i];
        }
        else
        {
          val = m_points[i].value;
        }
      }
      scalars->SetValue( i, val );
    }
    polydata->GetPointData()->SetScalars( scalars );
    polydata->Modified();
    scalars->Delete();
  }
}

// return index of the point if found, otherwise returns -1
int LayerPointSet::FindPoint( double* ras, double tolerance )
{
  double dt = tolerance;
  if ( dt < 0 )
  {
    double* voxel_size = m_layerRef->GetWorldVoxelSize();
    dt = GetProperty()->GetRadius() * qMin( voxel_size[0], qMin( voxel_size[1], voxel_size[2] ) );
    dt = dt * dt;
  }
  for ( int i = 0; i < m_points.size(); i++ )
  {
    if ( vtkMath::Distance2BetweenPoints( m_points[i].pt, ras ) < dt )
    {
      return i;
    }
  }
  return -1;
}

// returns index of the point
int LayerPointSet::AddPoint( double* ras_in, double value, bool bNotToVoxelCenter )
{
  int nRet;
  int dim[3];
  double ras[3], vs[3];
  m_layerRef->GetVolumeInfo(dim, vs);
  double min_tor2 = qMin(vs[0], qMin(vs[1], vs[2]))*3;
  min_tor2 *= min_tor2;
  if ( !bNotToVoxelCenter && GetProperty()->GetSnapToVoxelCenter() )
  {
    m_layerRef->SnapToVoxelCenter( ras_in, ras );
  }
  else
  {
    ras[0] = ras_in[0];
    ras[1] = ras_in[1];
    ras[2] = ras_in[2];
  }

  if ( m_points.size() < 2 || !GetProperty()->GetShowSpline())
  {
    ControlPoint p;
    p.pt[0] = ras[0];
    p.pt[1] = ras[1];
    p.pt[2] = ras[2];
    p.value = value;
    m_points.push_back( p );
    nRet = m_points.size() - 1;
  }
  else
  {
    // first find the closest point
    double dist = 1e20;
    int n = 0;
    for ( int i = 0; i < m_points.size(); i++ )
    {
      double temp = vtkMath::Distance2BetweenPoints( ras, m_points[i].pt );
      if ( temp < dist )
      {
        n = i;
        dist = temp;
      }
    }

    // then determine where to insert the point
    int n1, n2;
    if ( n == 0 )
    {
      n1 = 0;
      n2 = 1;
    }
    else
    {
      n1 = n-1;
      n2 = n;
    }

    double d1 = vtkMath::Distance2BetweenPoints( ras, m_points[n1].pt );
    double d2 = vtkMath::Distance2BetweenPoints( ras, m_points[n2].pt );
    double d3 = vtkMath::Distance2BetweenPoints( m_points[n1].pt, m_points[n2].pt );

    if ( d3 >= d1 && d3 >= d2)
    {
      n = n2;
    }
    else if ( d1 < d2 )
    {
      n = n1;
    }
    else
    {
      n = n2 + 1;
    }

    double d0 = vtkMath::Distance2BetweenPoints( ras, m_points[0].pt );
    double dn = vtkMath::Distance2BetweenPoints( ras, m_points[m_points.size()-1].pt );
//    if ((d0 > d1 && d0 < d2) || (d0 < d2 && d0 > d1) )
//      n = 0;
//    else if ((dn > d1 && dn < d2) || (dn < d2 && dn > d1) )
//      n = m_points.size();

    double pt[3];
    pt[0] = (m_points[n1].pt[0] + m_points[n2].pt[0])/2;
    pt[1] = (m_points[n1].pt[1] + m_points[n2].pt[1])/2;
    pt[2] = (m_points[n1].pt[2] + m_points[n2].pt[2])/2;
    double d4 = vtkMath::Distance2BetweenPoints( ras, pt );
    if (d4 > d3/3)
    {
      if ( d0 < dn )
        n = 0;
      else
        n = m_points.size();
    }

    ControlPoint p;
    p.pt[0] = ras[0];
    p.pt[1] = ras[1];
    p.pt[2] = ras[2];
    p.value = value;

    if ( n >= (int)m_points.size() )
    {
      m_points.push_back( p );
    }
    else
    {
      m_points.insert( m_points.begin() + n, p );
    }

    nRet = n;
  }

  SetModified();
  RebuildActors();
  emit PointAdded(nRet);
  return nRet;
}

bool LayerPointSet::RemovePoint( double* ras, double tolerance )
{
  Q_UNUSED(tolerance);
  return RemovePoint( FindPoint( ras ) );
}

bool LayerPointSet::RemovePoint( int nIndex )
{
  if ( nIndex < 0 || nIndex >= (int)m_points.size() )
  {
    return false;
  }

  m_points.erase( m_points.begin() + nIndex );

  SetModified();

  RebuildActors();

  emit PointRemoved(nIndex);

  return true;
}

void LayerPointSet::UpdatePoint( int nIndex, double* ras, bool rebuildActor )
{
  if (m_points.size() <= nIndex)
    return;

  if ( GetProperty()->GetSnapToVoxelCenter() )
  {
    m_layerRef->SnapToVoxelCenter( ras, m_points[nIndex].pt );
  }
  else
  {
    m_points[nIndex].pt[0] = ras[0];
    m_points[nIndex].pt[1] = ras[1];
    m_points[nIndex].pt[2] = ras[2];
  }

  SetModified();

  if ( rebuildActor )
  {
    RebuildActors();
  }
}

void LayerPointSet::UpdatePoint(int nIndex, const QString &key, const QVariant &value)
{
  QVariantMap map = m_points[nIndex].info;
  map[key] = value;
  m_points[nIndex].info = map;
  if (m_mapEnhancedData.isEmpty())
    m_mapEnhancedData["data_type"] = "fs_pointset";

  if (key == "legacy_stat")
    m_points[nIndex].value = value.toDouble();
  SetModified();
}

void LayerPointSet::UpdateSnapToVoxelCenter()
{
  for ( int i = 0; i < m_points.size(); i++ )
  {
    UpdatePoint(i, m_points[i].pt, false );
  }
  this->RebuildActors();
}

void LayerPointSet::UpdateOpacity()
{
  for ( int i = 0; i < 3; i++ )
  {
    m_actorSlice[i]->GetProperty()->SetOpacity( GetProperty()->GetOpacity() );
    m_actorSplineSlice[i]->GetProperty()->SetOpacity(GetProperty()->GetOpacity() );
  }
  emit ActorUpdated();
}

void LayerPointSet::UpdateColorMap()
{
  for ( int i = 0; i < 3; i++ )
  {
    m_actorSlice[i]->GetProperty()->SetColor( GetProperty()->GetColor() );
    m_actorSplineSlice[i]->GetProperty()->SetColor( GetProperty()->GetSplineColor() );
  }

  m_actorBalls->GetProperty()->SetColor( GetProperty()->GetColor() );

  switch ( GetProperty()->GetColorMap() )
  {
  case LayerPropertyPointSet::SolidColor:
    m_mapper->ScalarVisibilityOff();
    m_actorSpline->GetProperty()->SetColor( GetProperty()->GetSplineColor() );
    if (m_actorSplineSlice[0]->GetMapper())
    {
      for (int i = 0; i < 3; i++)
      {
        m_actorSplineSlice[i]->GetMapper()->ScalarVisibilityOff();
        m_actorSplineSlice[i]->GetProperty()->SetColor(GetProperty()->GetSplineColor());
      }
    }
    break;
  case LayerPropertyPointSet::HeatScale:
    m_mapper->ScalarVisibilityOn();
    m_mapper->SetLookupTable( GetProperty()->GetHeatScaleLUT() );
    if (m_actorSplineSlice[0]->GetMapper())
    {
      for (int i = 0; i < 3; i++)
      {
        m_actorSplineSlice[i]->GetMapper()->ScalarVisibilityOn();
        m_actorSplineSlice[i]->GetMapper()->SetLookupTable( GetProperty()->GetHeatScaleLUT() );
      }
    }
    break;
  }
  emit ActorUpdated();
}

bool LayerPointSet::Rotate( std::vector<RotationElement>& rotations, wxWindow* wnd, wxCommandEvent& event )
{
  Q_UNUSED(rotations);
  Q_UNUSED(wnd);
  Q_UNUSED(event);
  m_points.clear();
  m_pointSetSource->LabelToPointSet( m_points, m_layerRef->GetSourceVolume() );
  RebuildActors();

  return true;
}

std::vector<double> LayerPointSet::GetPoints()
{
  std::vector<double> list;
  foreach (ControlPoint wp, m_points)
  {
    list.push_back(wp.pt[0]);
    list.push_back(wp.pt[1]);
    list.push_back(wp.pt[2]);
  }

  return list;
}

void LayerPointSet::GetPoint(int nIndex, double *pt_out)
{
  if (nIndex < m_points.size())
  {
    pt_out[0] = m_points[nIndex].pt[0];
    pt_out[1] = m_points[nIndex].pt[1];
    pt_out[2] = m_points[nIndex].pt[2];
  }
}

ControlPoint LayerPointSet::GetPoint(int nIndex)
{
  if (nIndex < m_points.size())
    return m_points[nIndex];
  else
    return ControlPoint();
}

int LayerPointSet::GetNumberOfPoints()
{
  return m_points.size();
}

bool LayerPointSet::GetCentroidPosition(double *pos)
{
  m_pointSetSource->UpdateLabel( m_points, m_layerRef->GetSourceVolume() );
  if (m_pointSetSource->GetCentroidRASPosition(pos, m_layerRef->GetSourceVolume()))
  {
    m_layerRef->RASToTarget(pos, pos);
    return true;
  }
  else
    return false;
}

bool LayerPointSet::IsEnhanced()
{
  return !m_mapEnhancedData.isEmpty();
}

double LayerPointSet::GetEndPointDistance()
{
  double val = 0;
  if (m_points.size() > 1)
  {
    val = MyUtils::GetDistance<double>(m_points.first().pt, m_points.last().pt);
  }
  return val;
}

vtkPoints* LayerPointSet::GetSplinedPoints()
{
  return m_splinedPoints;
}

void LayerPointSet::GetNormalAtPoint(int nIndex, double *vnorm, int nPlane)
{
  double pt[3], pt0[3], pt1[3];
  GetPoint(nIndex, pt);
  int n = 0, n0, n1;
  double* vs = m_layerRef->GetWorldVoxelSize();
  double dTor2 = qMin(vs[0], qMin(vs[1], vs[2]))*0.02;
  dTor2 *= dTor2;
  double dMinDist2 = 1e8;
  for (int i = 0; i < m_splinedPoints->GetNumberOfPoints(); i++)
  {
    double dval = vtkMath::Distance2BetweenPoints(pt, m_splinedPoints->GetPoint(i));
    if (dval < dTor2)
    {
      n = i;
      break;
    }
    else if (dval < dMinDist2)
    {
      n = i;
      dMinDist2 = dval;
    }
  }
  if (n == 0)
  {
    n0 = 0;
    n1 = 1;
  }
  else if (n == m_splinedPoints->GetNumberOfPoints()-1)
  {
    n0 = m_splinedPoints->GetNumberOfPoints()-1;
    n1 = n0-1;
  }
  else
  {
    n0 = n-1;
    n1 = n+1;
  }
  m_splinedPoints->GetPoint(n0, pt0);
  m_splinedPoints->GetPoint(n1, pt1);
  double v[3], v2[3] = {0, 0, 0};
  for (int i = 0; i < 3; i++)
    v[i] = pt1[i] - pt0[i];
  vtkMath::Normalize(v);
  v2[nPlane] = 1;
  vtkMath::Cross(v, v2, vnorm);
}

void LayerPointSet::SetEnhancedData(const QString &key, const QVariant &val)
{
  m_mapEnhancedData[key] = val;
  SetModified();
}

QVariant LayerPointSet::GetEnhancedData(const QString &key)
{
  QVariant v;
  if (m_mapEnhancedData.contains(key))
    v = m_mapEnhancedData[key];
  return v;
}
