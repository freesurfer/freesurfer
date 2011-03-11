/**
 * @file  LayerPropertiesWayPoints.cxx
 * @brief Implementation for ROI layer properties.
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
#include <wx/wx.h>
#include <wx/ffile.h>
#include <wx/filename.h>
#include "LayerPropertiesWayPoints.h"
#include "vtkRGBAColorTransferFunction.h"
#include "LayerMRI.h"
#include "FSVolume.h"
#include "MyUtils.h"

using namespace std;

LayerPropertiesWayPoints::LayerPropertiesWayPoints () :
    LayerProperties()
{
  mOpacity = 0.7;
  mRGB[0] = 1;
  mRGB[1] = 0.1;
  mRGB[2] = 0.1;

  mRGBSpline[0] = 1;
  mRGBSpline[1] = 1;
  mRGBSpline[2] = 0;

  m_dRadius = 1;
  m_dSplineRadius = 0.5;
  m_nColorMap = SolidColor;
  m_layerScalar = NULL;
  m_nScalarSet = -1;

  m_dHeatScaleMin = 0;
  m_dHeatScaleMid = 0.5;
  m_dHeatScaleMax = 1;
  m_dHeatScaleOffset = 0;

  m_nScalarType = ScalarLayer;
  m_nType = WayPoints;

  m_bShowSpline = true;
  m_bSnapToVoxelCenter = false;
  
  m_lutHeatScale = vtkSmartPointer<vtkRGBAColorTransferFunction>::New();

  this->ColorMapChanged();
}

LayerPropertiesWayPoints::~LayerPropertiesWayPoints ()
{
  for ( size_t i = 0; i < m_scalarSets.size(); i++ )
  {
    if ( m_scalarSets[i].dValue )
      delete[] m_scalarSets[i].dValue;
  }
}


void LayerPropertiesWayPoints::ColorMapChanged ()
{
  switch ( m_nColorMap )
  {
  case SolidColor:
    break;
  case HeatScale:
    m_lutHeatScale->RemoveAllPoints();
    m_lutHeatScale->AddRGBAPoint( -m_dHeatScaleMax + m_dHeatScaleOffset, 0, 1, 1, 1 );
    m_lutHeatScale->AddRGBAPoint( -m_dHeatScaleMid + m_dHeatScaleOffset, 0, 0, 1, 1 );
    m_lutHeatScale->AddRGBAPoint( -m_dHeatScaleMin + m_dHeatScaleOffset, 0, 0, 1, 0 );
    m_lutHeatScale->AddRGBAPoint(  m_dHeatScaleMin + m_dHeatScaleOffset, 1, 0, 0, 0 );
    m_lutHeatScale->AddRGBAPoint(  m_dHeatScaleMid + m_dHeatScaleOffset, 1, 0, 0, 1 );
    m_lutHeatScale->AddRGBAPoint(  m_dHeatScaleMax + m_dHeatScaleOffset, 1, 1, 0, 1 );
    /* mHeatScaleTable->AddRGBAPoint(  mHeatScaleMinThreshold, 0, 0, 0, 0 );
     mHeatScaleTable->AddRGBAPoint(  mHeatScaleMinThreshold+0.0001, 0, 0, 0, 1 );
     for ( int i = 0; i < n; i++ )
      mHeatScaleTable->AddRGBAPoint( mHeatScaleMinThreshold +0.0001+ i*(mHeatScaleMidThreshold-mHeatScaleMinThreshold)/n,
                1.0*i/n, 0, 0, 1 );
     mHeatScaleTable->AddRGBAPoint(  mHeatScaleMidThreshold, 1, 0, 0, 1 );
     for ( int i = 0; i < n; i++ )
      mHeatScaleTable->AddRGBAPoint( mHeatScaleMidThreshold + i*(mHeatScaleMaxThreshold-mHeatScaleMidThreshold)/n,
                1, 1.0*i/n, 0, 1 );
     mHeatScaleTable->AddRGBAPoint(  mHeatScaleMaxThreshold, 1, 1, 0, 1 );*/
    m_lutHeatScale->Build();
    break;
  }

  // Notify the layers that use the color map stuff.
  this->SendBroadcast( "ColorMapChanged", NULL );
}

vtkRGBAColorTransferFunction* LayerPropertiesWayPoints::GetHeatScaleLUT()
{
  return m_lutHeatScale;
}

void LayerPropertiesWayPoints::SetColor ( double r, double g, double b )
{
  mRGB[0] = r;
  mRGB[1] = g;
  mRGB[2] = b;
  this->ColorMapChanged();
}

void LayerPropertiesWayPoints::SetSplineColor ( double r, double g, double b )
{
  mRGBSpline[0] = r;
  mRGBSpline[1] = g;
  mRGBSpline[2] = b;
  this->ColorMapChanged();
}

double LayerPropertiesWayPoints::GetOpacity() const
{
  return mOpacity;
}

void LayerPropertiesWayPoints::SetOpacity( double opacity )
{
  if ( mOpacity != opacity )
  {
    mOpacity = opacity;
    if ( mOpacity >= 1 )
      mOpacity = 0.99999;   // strange vtk behavior. If set to 1, it will be transparent!
    this->SendBroadcast( "OpacityChanged", this );
  }
}

void LayerPropertiesWayPoints::SetRadius( double r )
{
  if ( m_dRadius != r )
  {
    m_dRadius = r;
    this->SendBroadcast( "RadiusChanged", this );
  }
}

void LayerPropertiesWayPoints::SetSplineRadius( double r )
{
  if ( m_dSplineRadius != r )
  {
    m_dSplineRadius = r;
    this->SendBroadcast( "RadiusChanged", this );
  }
}

void LayerPropertiesWayPoints::SetColorMap( int map )
{
  if ( m_nColorMap != map )
  {
    m_nColorMap = map;
    this->ColorMapChanged();
  }
}

void LayerPropertiesWayPoints::SetHeatScaleMin( double dValue )
{
  if ( m_dHeatScaleMin != dValue )
  {
    m_dHeatScaleMin = dValue;
    this->ColorMapChanged();
  }
}

void LayerPropertiesWayPoints::SetHeatScaleMid( double dValue )
{
  if ( m_dHeatScaleMid != dValue )
  {
    m_dHeatScaleMid = dValue;
    this->ColorMapChanged();
  }
}

void LayerPropertiesWayPoints::SetHeatScaleMax( double dValue )
{
  if ( m_dHeatScaleMax != dValue )
  {
    m_dHeatScaleMax = dValue;
    this->ColorMapChanged();
  }
}

void LayerPropertiesWayPoints::SetHeatScaleOffset( double dValue )
{
  if ( m_dHeatScaleOffset != dValue )
  {
    m_dHeatScaleOffset = dValue;
    this->ColorMapChanged();
  }
}

void LayerPropertiesWayPoints::SetScalarLayer( LayerMRI* layer )
{
  if ( m_nScalarType != ScalarLayer || m_layerScalar != layer )
  {
    m_nScalarType = ScalarLayer;
    m_layerScalar = layer;
    if ( layer )
      layer->AddListener( this );
    UpdateScalarValues();
    this->ColorMapChanged();
    this->SendBroadcast( "ScalarLayerChanged", this );
  }
}

void LayerPropertiesWayPoints::SetScalarSet( int n )
{
  if ( n < (int)m_scalarSets.size() && n >= 0 )
  {
    m_nScalarType = ScalarSet;
    m_nScalarSet = n;
    UpdateScalarValues();
    this->ColorMapChanged();
    this->SendBroadcast( "ScalarSetChanged", this );
  }
}

void LayerPropertiesWayPoints::UpdateScalarValues()
{
  double fMin = GetScalarMinValue();
  double fMax = GetScalarMaxValue();
  m_dHeatScaleMin = fMin;
  m_dHeatScaleMid = ( fMax - fMin ) / 2;
  m_dHeatScaleMax = fMax;
  m_dHeatScaleOffset = 0;
}

double LayerPropertiesWayPoints::GetScalarMinValue()
{
  double min = 0;
  if ( (m_nScalarType == ScalarLayer) && m_layerScalar )
    min = m_layerScalar->GetSourceVolume()->GetMinValue();
  else if ( m_nScalarType == ScalarSet )
    min = GetActiveScalarSet().dMin;

  return min;
}

double LayerPropertiesWayPoints::GetScalarMaxValue()
{
  double max = 1;
  if ( (m_nScalarType == ScalarLayer) && m_layerScalar )
    max = m_layerScalar->GetSourceVolume()->GetMaxValue();
  else if ( m_nScalarType == ScalarSet )
    max = GetActiveScalarSet().dMax;
  
  return max;
}

void LayerPropertiesWayPoints::DoListenToMessage ( std::string const iMessage, 
						   void* iData, void* sender )
{
  if ( iMessage == "LayerObjectDeleted" && iData == m_layerScalar )
  {
    m_layerScalar = NULL;
  }
}

bool LayerPropertiesWayPoints::LoadScalarsFromFile( const char* filename )
{
  wxFFile file( wxString::FromAscii(filename) );
  wxString strg;
  if ( !file.ReadAll( &strg ) )
    return false;
  strg.Replace( _("\n"), _(" ") );
  wxArrayString ar = MyUtils::SplitString( strg, _(" ") );

  std::vector<double> values;
  double val;
  for ( size_t i = 0; i < ar.GetCount(); i++ )
  {
    ar[i].ToDouble( &val );
    values.push_back( val );
  }
  ScalarValues sv;
  sv.nNum = values.size();
  sv.dValue = new double[sv.nNum];
  sv.dMin = 1e8;
  sv.dMax = -1e8;
  for ( int i = 0; i < sv.nNum; i++ )
  {
    sv.dValue[i] = values[i];
    if ( sv.dValue[i] < sv.dMin )
      sv.dMin = sv.dValue[i];
    else if ( sv.dValue[i] > sv.dMax )
      sv.dMax = sv.dValue[i];
  }
  sv.strName = filename; // wxFileName( filename ).GetName().c_str();

  if ( sv.nNum < 2 )
    return false;

  m_scalarSets.push_back( sv );
  SetScalarSet( m_scalarSets.size() - 1 );

  return true;
}

void LayerPropertiesWayPoints::SetSnapToVoxelCenter( bool bSnag )
{
  m_bSnapToVoxelCenter = bSnag;
  
  this->SendBroadcast( "SnapToVoxelCenterChanged", this );
}

void LayerPropertiesWayPoints::SetShowSpline( bool bSpline )
{
  m_bShowSpline = bSpline;
  
  this->SendBroadcast( "SplineVisibilityChanged", this );
}

void LayerPropertiesWayPoints::SetType( int nType )
{
  m_nType = nType;
}
