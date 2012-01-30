/**
 * @file  LayerPropertyPointSet.cxx
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
 *    $Author: rpwang $
 *    $Date: 2012/01/30 20:57:05 $
 *    $Revision: 1.8 $
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


#include "LayerPropertyPointSet.h"
#include "vtkRGBAColorTransferFunction.h"
#include "LayerMRI.h"
#include "FSVolume.h"
#include "MyUtils.h"
#include <QFile>
#include <QStringList>
#include <QTextStream>

//using namespace std;

LayerPropertyPointSet::LayerPropertyPointSet (QObject* parent) :
  LayerProperty(parent)
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
  m_nType = WayPoint;

  m_bShowSpline = true;
  m_bSnapToVoxelCenter = false;

  m_lutHeatScale = vtkSmartPointer<vtkRGBAColorTransferFunction>::New();

  connect(this, SIGNAL(SnapToVoxelCenterChanged(bool)), this, SIGNAL(PropertyChanged()));
  connect(this, SIGNAL(SplineVisibilityChanged(bool)), this, SIGNAL(PropertyChanged()));
  connect(this, SIGNAL(ColorMapChanged()), this, SIGNAL(PropertyChanged()));
  connect(this, SIGNAL(OpacityChanged(double)), this, SIGNAL(PropertyChanged()));
  connect(this, SIGNAL(RadiusChanged(double)), this, SIGNAL(PropertyChanged()));
  connect(this, SIGNAL(SplineRadiusChanged(double)), this, SIGNAL(PropertyChanged()));
}

LayerPropertyPointSet::~LayerPropertyPointSet ()
{
  for ( size_t i = 0; i < m_scalarSets.size(); i++ )
  {
    if ( m_scalarSets[i].dValue )
    {
      delete[] m_scalarSets[i].dValue;
    }
  }
}

void LayerPropertyPointSet::SetColorMapChanged ()
{
  switch ( m_nColorMap )
  {
  case SolidColor:
    break;
  case HeatScale:
    m_lutHeatScale->RemoveAllPoints();
    /*
    m_lutHeatScale->AddRGBAPoint( -m_dHeatScaleMax + m_dHeatScaleOffset, 0, 1, 1, 1 );
    m_lutHeatScale->AddRGBAPoint( -m_dHeatScaleMid + m_dHeatScaleOffset, 0, 0, 1, 1 );
    m_lutHeatScale->AddRGBAPoint( -m_dHeatScaleMin + m_dHeatScaleOffset, 0, 0, 1, 0 );
    */
    m_lutHeatScale->AddRGBAPoint(  m_dHeatScaleMin + m_dHeatScaleOffset, 1, 0, 0, 1 );
    m_lutHeatScale->AddRGBAPoint(  m_dHeatScaleMid + m_dHeatScaleOffset, 1, 0.5, 0, 1 );
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
  emit ColorMapChanged();
}

vtkRGBAColorTransferFunction* LayerPropertyPointSet::GetHeatScaleLUT()
{
  return m_lutHeatScale;
}

void LayerPropertyPointSet::SetColor ( double r, double g, double b )
{
  mRGB[0] = r;
  mRGB[1] = g;
  mRGB[2] = b;
  this->SetColorMapChanged();
}

void LayerPropertyPointSet::SetSplineColor ( double r, double g, double b )
{
  mRGBSpline[0] = r;
  mRGBSpline[1] = g;
  mRGBSpline[2] = b;
  this->SetColorMapChanged();
}

double LayerPropertyPointSet::GetOpacity() const
{
  return mOpacity;
}

void LayerPropertyPointSet::SetOpacity( double opacity )
{
  if ( mOpacity != opacity )
  {
    mOpacity = opacity;
    if ( mOpacity >= 1 )
    {
      mOpacity = 0.99999;  // strange vtk behavior. If set to 1, it will be transparent!
    }
    emit OpacityChanged( opacity );
  }
}

void LayerPropertyPointSet::SetRadius( double r )
{
  if ( m_dRadius != r )
  {
    m_dRadius = r;
    emit RadiusChanged( r );
  }
}

void LayerPropertyPointSet::SetSplineRadius( double r )
{
  if ( m_dSplineRadius != r )
  {
    m_dSplineRadius = r;
    emit SplineRadiusChanged( r );
  }
}

void LayerPropertyPointSet::SetColorMap( int map )
{
  if ( m_nColorMap != map )
  {
    m_nColorMap = map;
    this->SetColorMapChanged();
  }
}

void LayerPropertyPointSet::SetHeatScaleMin( double dValue )
{
  if ( m_dHeatScaleMin != dValue )
  {
    m_dHeatScaleMin = dValue;
    this->SetColorMapChanged();
  }
}

void LayerPropertyPointSet::SetHeatScaleMid( double dValue )
{
  if ( m_dHeatScaleMid != dValue )
  {
    m_dHeatScaleMid = dValue;
    this->SetColorMapChanged();
  }
}

void LayerPropertyPointSet::SetHeatScaleMax( double dValue )
{
  if ( m_dHeatScaleMax != dValue )
  {
    m_dHeatScaleMax = dValue;
    this->SetColorMapChanged();
  }
}

void LayerPropertyPointSet::SetHeatScaleOffset( double dValue )
{
  if ( m_dHeatScaleOffset != dValue )
  {
    m_dHeatScaleOffset = dValue;
    this->SetColorMapChanged();
  }
}

void LayerPropertyPointSet::SetScalarLayer( LayerMRI* layer )
{
  if ( m_nScalarType != ScalarLayer || m_layerScalar != layer )
  {
    m_nScalarType = ScalarLayer;
    m_layerScalar = layer;
//   if ( layer )
//     layer->AddListener( this );
    UpdateScalarValues();
    this->SetColorMapChanged();
    emit ScalarLayerChanged( layer );
  }
}

void LayerPropertyPointSet::SetScalarSet( int n )
{
  if ( n < (int)m_scalarSets.size() && n >= 0 )
  {
    m_nScalarType = ScalarSet;
    m_nScalarSet = n;
    UpdateScalarValues();
    this->SetColorMapChanged();
    emit ScalarSetChanged();
  }
}

void LayerPropertyPointSet::UpdateScalarValues()
{
  double fMin = GetScalarMinValue();
  double fMax = GetScalarMaxValue();
  m_dHeatScaleMin = fMin;
  m_dHeatScaleMid = ( fMax - fMin ) / 2;
  m_dHeatScaleMax = fMax;
  m_dHeatScaleOffset = 0;
}

double LayerPropertyPointSet::GetScalarMinValue()
{
  double min = 0;
  if ( (m_nScalarType == ScalarLayer) && m_layerScalar )
  {
    min = m_layerScalar->GetSourceVolume()->GetMinValue();
  }
  else if ( m_nScalarType == ScalarSet )
  {
    min = GetActiveScalarSet().dMin;
  }

  return min;
}

double LayerPropertyPointSet::GetScalarMaxValue()
{
  double max = 1;
  if ( (m_nScalarType == ScalarLayer) && m_layerScalar )
  {
    max = m_layerScalar->GetSourceVolume()->GetMaxValue();
  }
  else if ( m_nScalarType == ScalarSet )
  {
    max = GetActiveScalarSet().dMax;
  }

  return max;
}

bool LayerPropertyPointSet::LoadScalarsFromFile( const QString& filename )
{
  QFile file(filename);
  if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
  {
    cerr << qPrintable(file.errorString()) << "\n";;
    return false;
  }

  QTextStream in(&file);
  std::vector<double> values;
  while (!in.atEnd())
  {
    QStringList strgs = in.readLine().split( " ", QString::SkipEmptyParts );
    for ( int i = 0; i < strgs.size(); i++ )
    {
      values.push_back( strgs[i].toDouble() );
    }
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
    {
      sv.dMin = sv.dValue[i];
    }
    else if ( sv.dValue[i] > sv.dMax )
    {
      sv.dMax = sv.dValue[i];
    }
  }
  sv.strName = filename;

  if ( sv.nNum < 2 )
  {
    return false;
  }

  m_scalarSets.push_back( sv );
  SetScalarSet( m_scalarSets.size() - 1 );

  return true;
}

void LayerPropertyPointSet::SetSnapToVoxelCenter( bool bSnag )
{
  m_bSnapToVoxelCenter = bSnag;

  emit SnapToVoxelCenterChanged( bSnag );
}

void LayerPropertyPointSet::SetShowSpline( bool bSpline )
{
  m_bShowSpline = bSpline;

  emit SplineVisibilityChanged( bSpline );
}

void LayerPropertyPointSet::SetType( int nType )
{
  m_nType = nType;
}
