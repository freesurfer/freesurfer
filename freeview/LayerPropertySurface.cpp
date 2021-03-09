/**
 * @brief Implementation for surface layer properties.
 *
 * In 2D, the MRI is viewed as a single slice, and controls are
 * provided to change the color table and other viewing options. In
 * 3D, the MRI is viewed in three planes in 3D space, with controls to
 * move each plane axially.
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


#include <assert.h>
#include "LayerPropertySurface.h"
#include "vtkLookupTable.h"
#include "vtkRGBAColorTransferFunction.h"
#include "vtkMath.h"
#include "FSSurface.h"
#include <QDebug>

using namespace std;

LayerPropertySurface::LayerPropertySurface ( QObject* parent ) :
  LayerProperty( parent ),
  m_dOpacity( 1 ),
  m_nEdgeThickness( 2 ),
  m_nVectorPointSize( 3 ),
  m_dThresholdMidPoint( 0 ),
  m_dThresholdSlope( 10 ),
  m_nCurvatureMap( CM_Threshold ),
  m_nSurfaceRenderMode( SM_Surface ),
  m_bShowVertices( false ),
  m_nVertexPointSize( 3 ),
  m_nMeshColorMap( 0 ),
  m_surface( NULL ),
  m_bShowOverlay(true),
  m_bShowAnnotation(true),
  m_bUseSurfaceColorOn2D(false),
  m_nZOrderLabel(3),
  m_nZOrderAnnotation(2),
  m_nZOrderOverlay(1)
{
  m_lutCurvature = vtkSmartPointer<vtkRGBAColorTransferFunction>::New();

  blockSignals( true );
  SetBinaryColor( 0.4, 0.4, 0.4 );
  double c1[3] = { 0, 1, 0 }, c2[3] = { 1, 0, 0 };
  SetThresholdColor( c1, c2 );
  SetEdgeColor( 1, 1, 0 );
  SetVectorColor( 1, 0.75, 0 );
  SetVertexColor( 0., 1.0, 1.0 );
  SetMeshColor( 0.75, 0.75, 0.75 );
  blockSignals( false );
  for ( int i = 0; i < 3; i++ )
  {
    m_dPosition[i] = 0;
  }

  connect( this, SIGNAL(ColorMapChanged()), this, SIGNAL(PropertyChanged()) );
  connect( this, SIGNAL(DisplayModeChanged()), this, SIGNAL(PropertyChanged()) );
  connect( this, SIGNAL(EdgeThicknessChanged(int)), this, SIGNAL(PropertyChanged()) );
  connect( this, SIGNAL(MeshRenderChanged()), this, SIGNAL(PropertyChanged()) );
  connect( this, SIGNAL(OpacityChanged(double)), this, SIGNAL(PropertyChanged()) );
  connect( this, SIGNAL(PositionChanged()), this, SIGNAL(PropertyChanged()) );
  connect( this, SIGNAL(RenderModeChanged(int)), this, SIGNAL(PropertyChanged()) );
  connect( this, SIGNAL(VectorPointSizeChanged(int)), this, SIGNAL(PropertyChanged()) );
  connect( this, SIGNAL(VertexRenderChanged()), this, SIGNAL(PropertyChanged()) );
}

LayerPropertySurface::~LayerPropertySurface ()
{}

void LayerPropertySurface::SetColorMapChanged()
{
  BuildCurvatureLUT( m_lutCurvature, m_nCurvatureMap );
  // Notify the layers that use the color map stuff.
  emit ColorMapChanged();
}

void LayerPropertySurface::RebuildCurvatureLUT()
{
  SetColorMapChanged();
}

QVariantMap LayerPropertySurface::GetFullSettings()
{
  QVariantMap map;
  map["Opacity"]  = m_dOpacity;
  map["RGB_R"]    = m_dRGB[0];
  map["RGB_G"]    = m_dRGB[1];
  map["RGB_B"]    = m_dRGB[2];
  map["RGBThresholdHigh_R"]   = m_dRGBThresholdHigh[0];
  map["RGBThresholdHigh_G"]   = m_dRGBThresholdHigh[1];
  map["RGBThresholdHigh_B"]   = m_dRGBThresholdHigh[2];
  map["RGBThresholdLow_R"]    = m_dRGBThresholdLow[0];
  map["RGBThresholdLow_G"]    = m_dRGBThresholdLow[1];
  map["RGBThresholdLow_B"]    = m_dRGBThresholdLow[2];
  map["RGBEdge_R"]  = m_dRGBEdge[0];
  map["RGBEdge_G"]  = m_dRGBEdge[1];
  map["RGBEdge_B"]  = m_dRGBEdge[2];
  map["RGBVector_R"]  = m_dRGBVector[0];
  map["RGBVector_G"]  = m_dRGBVector[1];
  map["RGBVector_B"]  = m_dRGBVector[2];
  map["RGBMesh_R"]  = m_dRGBMesh[0];
  map["RGBMesh_G"]  = m_dRGBMesh[1];
  map["RGBMesh_B"]  = m_dRGBMesh[2];
  map["EdgeThickness"] = m_nEdgeThickness;
  map["VectorPointSize"] = m_nVectorPointSize;

  map["ThresholdMidPoint"]  = m_dThresholdMidPoint;
  map["ThresholdSlope"]     =  m_dThresholdSlope;

  map["Position_X"] = m_dPosition[0];
  map["Position_Y"] = m_dPosition[1];
  map["Position_Z"] = m_dPosition[2];

  map["CurvatureMap"] = m_nCurvatureMap;

  map["SurfaceRenderMode"] = m_nSurfaceRenderMode;

  map["ShowVertices"] = m_bShowVertices;
  map["RGBVertex_R"] = m_dRGBVertex[0];
  map["RGBVertex_G"] = m_dRGBVertex[1];
  map["RGBVertex_B"] = m_dRGBVertex[2];
  map["VertexPointSize"] = m_nVertexPointSize;

  map["MeshColorMap"] = m_nMeshColorMap;

  return map;
}

void LayerPropertySurface::RestoreFullSettings(const QVariantMap &map)
{
  if (map.contains("Opacity"))
  {
    SetOpacity(map["Opacity"].toDouble());
  }
  if (map.contains("RGB_R"))
  {
    m_dRGB[0] = map["RGB_R"].toDouble();
    m_dRGB[1] = map["RGB_G"].toDouble();
    m_dRGB[2] = map["RGB_B"].toDouble();
  }
  if (map.contains("RGBThresholdHigh_R"))
  {
    m_dRGBThresholdHigh[0] = map["RGBThresholdHigh_R"].toDouble();
    m_dRGBThresholdHigh[1] = map["RGBThresholdHigh_G"].toDouble();
    m_dRGBThresholdHigh[2] = map["RGBThresholdHigh_B"].toDouble();
  }
  if (map.contains("RGBThresholdLow_R"))
  {
    m_dRGBThresholdLow[0] = map["RGBThresholdLow_R"].toDouble();
    m_dRGBThresholdLow[1] = map["RGBThresholdLow_G"].toDouble();
    m_dRGBThresholdLow[2] = map["RGBThresholdLow_B"].toDouble();
  }
  if (map.contains("RGBEdge_R"))
  {
    m_dRGBEdge[0] = map["RGBEdge_R"].toDouble();
    m_dRGBEdge[1] = map["RGBEdge_G"].toDouble();
    m_dRGBEdge[2] = map["RGBEdge_B"].toDouble();
  }
  if (map.contains("RGBVector_R"))
  {
    m_dRGBVector[0] = map["RGBVector_R"].toDouble();
    m_dRGBVector[1] = map["RGBVector_G"].toDouble();
    m_dRGBVector[2] = map["RGBVector_B"].toDouble();
  }
  if (map.contains("RGBMesh_R"))
  {
    m_dRGBMesh[0] = map["RGBMesh_R"].toDouble();
    m_dRGBMesh[1] = map["RGBMesh_G"].toDouble();
    m_dRGBMesh[2] = map["RGBMesh_B"].toDouble();
  }
  if (map.contains("RGBVertex_R"))
  {
    m_dRGBVertex[0] = map["RGBVertex_R"].toDouble();
    m_dRGBVertex[1] = map["RGBVertex_G"].toDouble();
    m_dRGBVertex[2] = map["RGBVertex_B"].toDouble();
  }
  if (map.contains("EdgeThickness"))
    SetEdgeThickness(map["EdgeThickness"].toInt());
  if (map.contains("VectorPointSize"))
    SetVectorPointSize(map["VectorPointSize"].toInt());
  if (map.contains("ThresholdMidPoint"))
    m_dThresholdMidPoint = map["ThresholdMidPoint"].toDouble();
  if (map.contains("ThresholdSlope"))
    m_dThresholdSlope = map["ThresholdSlope"].toDouble();

  if (map.contains("CurvatureMap"))
    m_nCurvatureMap = map["CurvatureMap"].toInt();
  if (map.contains("SurfaceRenderMode"))
    SetSurfaceRenderMode(map["SurfaceRenderMode"].toInt());

  if (map.contains("Position_X"))
  {
    double pos[3] = {map["Position_X"].toDouble(),
                     map["Position_Y"].toDouble(),
                     map["Position_Z"].toDouble()};
    SetPosition( pos );
  }
  if (map.contains("ShowVertices"))
    ShowVertices(map["ShowVertices"].toBool());

  if (map.contains("VertexPointSize"))
    SetVertexPointSize(map["VertexPointSize"].toInt());
  if (map.contains("MeshColorMap"))
    SetMeshColorMap(map["MeshColorMap"].toInt());

  SetColorMapChanged();
}

void LayerPropertySurface::BuildCurvatureLUT( vtkRGBAColorTransferFunction* lut, int nMap )
{
  double hiRGB[3];
  vtkMath::RGBToHSV( m_dRGB, hiRGB );
  hiRGB[2] *= 1.5;
  if ( hiRGB[2] > 1 )
  {
    hiRGB[2] /= (1.5*1.5);
  }
  vtkMath::HSVToRGB( hiRGB, hiRGB );
  lut->RemoveAllPoints();
  double dSlope = ( m_dThresholdSlope == 0 ? 1e8 : ( 1.0 / m_dThresholdSlope ) );
  switch ( nMap )
  {
  case CM_Threshold:
    lut->AddRGBAPoint( -dSlope + m_dThresholdMidPoint,
                       m_dRGBThresholdLow[0],
        m_dRGBThresholdLow[1],
        m_dRGBThresholdLow[2],
        1 );
    lut->AddRGBAPoint( m_dThresholdMidPoint, m_dRGB[0], m_dRGB[1], m_dRGB[2], 1 );
    lut->AddRGBAPoint( dSlope + m_dThresholdMidPoint,
                       m_dRGBThresholdHigh[0],
        m_dRGBThresholdHigh[1],
        m_dRGBThresholdHigh[2],
        1 );
    break;
  case CM_Binary:
    lut->AddRGBAPoint( m_dThresholdMidPoint,    hiRGB[0], hiRGB[1], hiRGB[2], 1 );
    lut->AddRGBAPoint( m_dThresholdMidPoint + 1e-8,    m_dRGB[0], m_dRGB[1], m_dRGB[2], 1 );
    break;
  default:
    lut->AddRGBAPoint( 0, m_dRGB[0], m_dRGB[1], m_dRGB[2], 1 );
    break;
  }

  lut->Build();
}

void LayerPropertySurface::SetSurfaceSource( FSSurface* surf )
{
  m_surface = surf;

  this->SetColorMapChanged();
}

vtkRGBAColorTransferFunction* LayerPropertySurface::GetCurvatureLUT() const
{
  return m_lutCurvature;
}

void LayerPropertySurface::SetCurvatureMap( int nMap )
{
  if ( m_nCurvatureMap != nMap )
  {
    m_nCurvatureMap = nMap;
    this->SetColorMapChanged();
  }
}

void LayerPropertySurface::SetThresholdMidPoint( double dValue )
{
  if ( m_dThresholdMidPoint != dValue )
  {
    m_dThresholdMidPoint = dValue;
    this->SetColorMapChanged();
  }
}

void LayerPropertySurface::SetThresholdSlope( double dValue )
{
  if ( m_dThresholdSlope != dValue )
  {
    m_dThresholdSlope = dValue;
    this->SetColorMapChanged();
  }
}


void LayerPropertySurface::SetThresholdColor( double* low, double* high )
{
  for ( int i = 0; i < 3; i++ )
  {
    m_dRGBThresholdLow[i] = low[i];
    m_dRGBThresholdHigh[i] = high[i];
  }
  this->SetColorMapChanged();
}

void LayerPropertySurface::SetBinaryColor ( double r, double g, double b )
{
  m_dRGB[0] = r;
  m_dRGB[1] = g;
  m_dRGB[2] = b;
  if ( m_surface )
  {
    if ( m_surface->IsCurvatureLoaded() )
    {}
    this->SetColorMapChanged();
  }
}

void LayerPropertySurface::SetEdgeColor ( double r, double g, double b )
{
  m_dRGBEdge[0] = r;
  m_dRGBEdge[1] = g;
  m_dRGBEdge[2] = b;
  emit ColorMapChanged();
  emit EdgeColorChanged();
}

void LayerPropertySurface::SetVectorColor ( double r, double g, double b )
{
  m_dRGBVector[0] = r;
  m_dRGBVector[1] = g;
  m_dRGBVector[2] = b;
  emit ColorMapChanged();
}

double LayerPropertySurface::GetOpacity() const
{
  return m_dOpacity;
}

void LayerPropertySurface::SetOpacity( double opacity )
{
  if ( m_dOpacity != opacity )
  {
    m_dOpacity = opacity;
    emit OpacityChanged( opacity );
  }
}

void LayerPropertySurface::SetEdgeThickness( int nThickness )
{
  if ( m_nEdgeThickness != nThickness )
  {
    m_nEdgeThickness = nThickness;
    emit EdgeThicknessChanged( nThickness );
  }
}

void LayerPropertySurface::SetVectorPointSize( int nSize )
{
  if ( m_nVectorPointSize != nSize )
  {
    m_nVectorPointSize = nSize;
    emit VectorPointSizeChanged( nSize );
  }
}

void LayerPropertySurface::SetSurfaceRenderMode( int nMode )
{
  if ( m_nSurfaceRenderMode != nMode )
  {
    m_nSurfaceRenderMode = nMode;
    emit RenderModeChanged( nMode );
  }
}

void LayerPropertySurface::ShowVertices( bool bShow )
{
  if ( m_bShowVertices != bShow )
  {
    m_bShowVertices = bShow;
    emit VertexRenderChanged();
  }
}

void LayerPropertySurface::SetVertexPointSize( int nSize )
{
  if ( m_nVertexPointSize != nSize )
  {
    m_nVertexPointSize = nSize;
    emit VertexRenderChanged();
  }
}

void LayerPropertySurface::SetVertexColor ( double r, double g, double b )
{
  m_dRGBVertex[0] = r;
  m_dRGBVertex[1] = g;
  m_dRGBVertex[2] = b;
  emit VertexRenderChanged();;
}

void LayerPropertySurface::SetMeshColor( double r, double g, double b )
{
  m_dRGBMesh[0] = r;
  m_dRGBMesh[1] = g;
  m_dRGBMesh[2] = b;
  emit MeshRenderChanged();
}

void LayerPropertySurface::SetMeshColorMap( int nMap )
{
  m_nMeshColorMap = nMap;
  emit MeshRenderChanged();
}

void LayerPropertySurface::SetPosition( double* p )
{
  double dp[3];
  for ( int i = 0; i < 3; i++ )
  {
    dp[i] = p[i] - m_dPosition[i];
    m_dPosition[i] = p[i];
  }

  emit PositionChanged();
  emit PositionChanged(dp[0], dp[1], dp[2]);
}


void LayerPropertySurface::SetBinaryColor( const QColor& c )
{
  SetBinaryColor( c.redF(), c.greenF(), c.blueF() );
}

void LayerPropertySurface::SetEdgeColor( const QColor& c )
{
  SetEdgeColor( c.redF(), c.greenF(), c.blueF() );
}

void LayerPropertySurface::SetVectorColor( const QColor& c )
{
  SetVectorColor( c.redF(), c.greenF(), c.blueF() );
}

void LayerPropertySurface::SetVertexColor( const QColor& c )
{
  SetVertexColor( c.redF(), c.greenF(), c.blueF() );
}

void LayerPropertySurface::SetMeshColor( const QColor& c )
{
  SetMeshColor( c.redF(), c.greenF(), c.blueF() );
}

void LayerPropertySurface::SetShowOverlay(bool bShow)
{
  if (bShow != m_bShowOverlay)
  {
    m_bShowOverlay = bShow;
    emit OverlayChanged();
  }
}

void LayerPropertySurface::SetShowAnnotation(bool bShow)
{
  if (bShow != m_bShowAnnotation)
  {
    m_bShowAnnotation = bShow;
    emit OverlayChanged();
  }
}

void LayerPropertySurface::SetUseSurfaceColorOn2D(bool bKeep)
{
  if (bKeep != m_bUseSurfaceColorOn2D)
  {
    m_bUseSurfaceColorOn2D = bKeep;
    emit ColorMapChanged();
  }
}

void LayerPropertySurface::SetZOrderAnnotation(int nOrder)
{
  if (m_nZOrderAnnotation != nOrder)
  {
    m_nZOrderAnnotation = nOrder;
    emit OverlayChanged();
  }
}

void LayerPropertySurface::SetZOrderLabel(int nOrder)
{
  if (m_nZOrderLabel != nOrder)
  {
    m_nZOrderLabel = nOrder;
    emit OverlayChanged();
  }
}

void LayerPropertySurface::SetZOrderOverlay(int nOrder)
{
  if (m_nZOrderOverlay != nOrder)
  {
    m_nZOrderOverlay = nOrder;
    emit OverlayChanged();
  }
}
