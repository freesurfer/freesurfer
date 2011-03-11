/**
 * @file  LayerPropertiesSurface.cxx
 * @brief Implementation for surface layer properties.
 *
 * In 2D, the MRI is viewed as a single slice, and controls are
 * provided to change the color table and other viewing options. In
 * 3D, the MRI is viewed in three planes in 3D space, with controls to
 * move each plane axially.
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


#include <assert.h>
#include "LayerPropertiesSurface.h"
#include "vtkLookupTable.h"
#include "vtkRGBAColorTransferFunction.h"
#include "vtkMath.h"
#include "FSSurface.h"
#include <wx/filename.h>

using namespace std;

LayerPropertiesSurface::LayerPropertiesSurface () :
    LayerProperties( ),
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
    m_surface( NULL )
{
  m_lutCurvature = vtkSmartPointer<vtkRGBAColorTransferFunction>::New();

  SetBinaryColor( 0.4, 0.4, 0.4 );
  double c1[3] = { 0, 1, 0 }, c2[3] = { 1, 0, 0 };
  SetThresholdColor( c1, c2 );
  SetEdgeColor( 1, 1, 0 );
  SetVectorColor( 1, 0.75, 0 );
  SetVertexColor( 0.75, 0.75, 0.75 );
  SetMeshColor( 0.75, 0.75, 0.75 );
  for ( int i = 0; i < 3; i++ )
    m_dPosition[i] = 0;
}

LayerPropertiesSurface::~LayerPropertiesSurface ()
{}

void LayerPropertiesSurface::ColorMapChanged()
{
  BuildCurvatureLUT( m_lutCurvature, m_nCurvatureMap );
  // Notify the layers that use the color map stuff.
  this->SendBroadcast( "ColorMapChanged", NULL );
}

void LayerPropertiesSurface::BuildCurvatureLUT( vtkRGBAColorTransferFunction* lut, int nMap )
{
  double hiRGB[3];
  vtkMath::RGBToHSV( m_dRGB, hiRGB );
  hiRGB[2] *= 1.5;
  if ( hiRGB[2] > 1 )
    hiRGB[2] /= (1.5*1.5);
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

void LayerPropertiesSurface::SetSurfaceSource( FSSurface* surf )
{
  m_surface = surf;

  this->ColorMapChanged();
}

vtkRGBAColorTransferFunction* LayerPropertiesSurface::GetCurvatureLUT() const
{
  return m_lutCurvature;
}

void LayerPropertiesSurface::SetCurvatureMap( int nMap )
{
  if ( m_nCurvatureMap != nMap )
  {
    m_nCurvatureMap = nMap;
    this->ColorMapChanged();
  }
}

void LayerPropertiesSurface::SetThresholdMidPoint( double dValue )
{
  if ( m_dThresholdMidPoint != dValue )
  {
    m_dThresholdMidPoint = dValue;
    this->ColorMapChanged();
  }
}

void LayerPropertiesSurface::SetThresholdSlope( double dValue )
{
  if ( m_dThresholdSlope != dValue )
  {
    m_dThresholdSlope = dValue;
    this->ColorMapChanged();
  }
}


void LayerPropertiesSurface::SetThresholdColor( double* low, double* high )
{
  for ( int i = 0; i < 3; i++ )
  {
    m_dRGBThresholdLow[i] = low[i];
    m_dRGBThresholdHigh[i] = high[i];
  }
  this->ColorMapChanged();
}

void LayerPropertiesSurface::SetBinaryColor ( double r, double g, double b )
{
  m_dRGB[0] = r;
  m_dRGB[1] = g;
  m_dRGB[2] = b;
  if ( m_surface )
  {
    if ( m_surface->IsCurvatureLoaded() )
    {}
    this->ColorMapChanged();
  }
}

void LayerPropertiesSurface::SetEdgeColor ( double r, double g, double b )
{
  m_dRGBEdge[0] = r;
  m_dRGBEdge[1] = g;
  m_dRGBEdge[2] = b;
  this->SendBroadcast( "ColorMapChanged", NULL );
}

void LayerPropertiesSurface::SetVectorColor ( double r, double g, double b )
{
  m_dRGBVector[0] = r;
  m_dRGBVector[1] = g;
  m_dRGBVector[2] = b;
  this->SendBroadcast( "ColorMapChanged", NULL );
}

double LayerPropertiesSurface::GetOpacity() const
{
  return m_dOpacity;
}

void LayerPropertiesSurface::SetOpacity( double opacity )
{
  if ( m_dOpacity != opacity )
  {
    m_dOpacity = opacity;
    this->SendBroadcast( "OpacityChanged", this );
  }
}

void LayerPropertiesSurface::SetEdgeThickness( int nThickness )
{
  if ( m_nEdgeThickness != nThickness )
  {
    m_nEdgeThickness = nThickness;
    this->SendBroadcast( "EdgeThicknessChanged", this );
  }
}

void LayerPropertiesSurface::SetVectorPointSize( int nSize )
{
  if ( m_nVectorPointSize != nSize )
  {
    m_nVectorPointSize = nSize;
    this->SendBroadcast( "VectorPointSizeChanged", this );
  }
}

void LayerPropertiesSurface::SetSurfaceRenderMode( int nMode )
{
  if ( m_nSurfaceRenderMode != nMode )
  {
    m_nSurfaceRenderMode = nMode;
    this->SendBroadcast( "SurfaceRenderModeChanged", this );
  }
}

void LayerPropertiesSurface::ShowVertices( bool bShow )
{
  if ( m_bShowVertices != bShow )
  {
    m_bShowVertices = bShow;
    this->SendBroadcast( "VertexRenderChanged", this );
  }
}

void LayerPropertiesSurface::SetVertexPointSize( int nSize )
{
  if ( m_nVertexPointSize != nSize )
  {
    m_nVertexPointSize = nSize;
    this->SendBroadcast( "VertexRenderChanged", this );
  }
}

void LayerPropertiesSurface::SetVertexColor ( double r, double g, double b )
{
  m_dRGBVertex[0] = r;
  m_dRGBVertex[1] = g;
  m_dRGBVertex[2] = b;
  this->SendBroadcast( "VertexRenderChanged", NULL );
}

void LayerPropertiesSurface::SetMeshColor( double r, double g, double b )
{
  m_dRGBMesh[0] = r;
  m_dRGBMesh[1] = g;
  m_dRGBMesh[2] = b;
  this->SendBroadcast( "MeshRenderChanged", NULL );
}

void LayerPropertiesSurface::SetMeshColorMap( int nMap )
{
  m_nMeshColorMap = nMap;
  this->SendBroadcast( "MeshRenderChanged", NULL );
}

void LayerPropertiesSurface::SetPosition( double* p )
{
  for ( int i = 0; i < 3; i++ )
    m_dPosition[i] = p[i];
  
  this->SendBroadcast( "PositionChanged", NULL );
}
