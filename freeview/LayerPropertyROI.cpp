/**
 * @brief Implementation for ROI layer properties.
 *
 * In 2D, the MRI is viewed as a single slice, and controls are
 * provided to change the color table and other viewing options. In
 * 3D, the MRI is viewed in three planes in 3D space, with controls to
 * move each plane axially.
 */
/*
 * Original Author: Kevin Teich
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
#include "LayerPropertyROI.h"
#include "vtkRGBAColorTransferFunction.h"
#include <QDebug>

using namespace std;

LayerPropertyROI::LayerPropertyROI ( QObject* parent) : LayerProperty( parent )
{
  mOpacity = 0.7;
  mRGB[0] = 1;
  mRGB[1] = 1;
  mRGB[2] = 0;
  m_dThreshold = 0;
  m_nColorCode = SolidColor;
  m_dValueRange[0] = 0;
  m_dValueRange[1] = 1;

  mLUTTable = vtkSmartPointer<vtkRGBAColorTransferFunction>::New();
  // mLUTTable->ClampingOff();
  this->SetColorMapChanged();

  connect( this, SIGNAL(ColorMapChanged()), this, SIGNAL(PropertyChanged()) );
  connect( this, SIGNAL(OpacityChanged(double)), this, SIGNAL(PropertyChanged()) );
  connect( this, SIGNAL(ThresholdChanged(double)), this, SIGNAL(PropertyChanged()));
}

LayerPropertyROI::~LayerPropertyROI ()
{}


vtkRGBAColorTransferFunction* LayerPropertyROI::GetLookupTable () const
{
  return mLUTTable;
}

void LayerPropertyROI::SetColorMapChanged()
{
  UpdateLUTTable();

  // Notify the layers that use the color map stuff.
  emit ColorMapChanged();
}

void LayerPropertyROI::UpdateLUTTable()
{
  assert( mLUTTable.GetPointer() );

  switch (m_nColorCode)
  {
  case SolidColor:
  {
    double range[2] = { m_dValueRange[0], m_dValueRange[1] };
    if (range[0] > 1)
      range[0] = 1;
    mLUTTable->RemoveAllPoints();
    mLUTTable->AddRGBAPoint( range[0]-0.0001, 0, 0, 0, 0 );
    mLUTTable->AddRGBAPoint( range[0],  mRGB[0], mRGB[1], mRGB[2], 1 );
    mLUTTable->AddRGBAPoint( range[0],  mRGB[0], mRGB[1], mRGB[2], 1 );
    mLUTTable->AddRGBAPoint( range[1],  mRGB[0], mRGB[1], mRGB[2], 1 );
    //      mLUTTable->AddRGBAPoint( range[1]+0.0001, 0, 0, 0, 0 );
  }

    //    mLUTTable->AddRGBAPoint( 1-0.0001, 0, 0, 0, 0 );
    //    mLUTTable->AddRGBAPoint( 1,  mRGB[0], mRGB[1], mRGB[2], 1 );
    //    mLUTTable->AddRGBAPoint( 1+0.0001, 0, 0, 0, 0 );

    break;
  case Heatscale:
    mLUTTable->RemoveAllPoints();
    mLUTTable->AddRGBAPoint( -m_dHeatscaleMax, 0, 1, 1, 1 );
    mLUTTable->AddRGBAPoint( -m_dHeatscaleMin, 0, 0, 1, 1 );
    mLUTTable->AddRGBAPoint(  0, 0, 0, 0, 0 );
    mLUTTable->AddRGBAPoint(  m_dHeatscaleMin, 1, 0, 0, 1 );
    mLUTTable->AddRGBAPoint(  m_dHeatscaleMax, 1, 1, 0, 1 );
    break;
  }
  mLUTTable->Build();
}

void LayerPropertyROI::SetColor ( double r, double g, double b )
{
  mRGB[0] = r;
  mRGB[1] = g;
  mRGB[2] = b;
  this->SetColorMapChanged();
}

void LayerPropertyROI::SetColorCode(int nCode)
{
  if (m_nColorCode != nCode)
  {
    m_nColorCode = nCode;
    SetColorMapChanged();
  }
}

void LayerPropertyROI::SetHeatscaleMax(double val)
{
  if (m_dHeatscaleMax != val)
  {
    m_dHeatscaleMax = val;
    SetColorMapChanged();
  }
}


void LayerPropertyROI::SetHeatscaleMin(double val)
{
  if (m_dHeatscaleMin != val)
  {
    m_dHeatscaleMin = val;
    SetColorMapChanged();
  }
}

void LayerPropertyROI::SetHeatscaleValues(double dMin, double dMax)
{
  m_dHeatscaleMin = dMin;
  m_dHeatscaleMax = dMax;
  SetColorMapChanged();
}

double LayerPropertyROI::GetOpacity() const
{
  return mOpacity;
}

void LayerPropertyROI::SetOpacity( double opacity )
{
  if ( mOpacity != opacity )
  {
    mOpacity = opacity;
    emit OpacityChanged( opacity );
  }
}

void LayerPropertyROI::SetThreshold(double th)
{
  if (th != m_dThreshold)
  {
    m_dThreshold = th;
    emit ThresholdChanged( th );
  }
}
