/**
 * @file  LayerPropertyROI.cxx
 * @brief Implementation for ROI layer properties.
 *
 * In 2D, the MRI is viewed as a single slice, and controls are
 * provided to change the color table and other viewing options. In
 * 3D, the MRI is viewed in three planes in 3D space, with controls to
 * move each plane axially.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2011/03/14 21:20:58 $
 *    $Revision: 1.3 $
 *
 * Copyright (C) 2007-2009,
 * The General Hospital Corporation (Boston, MA).
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 * Bug reports: analysis-bugs@nmr.mgh.harvard.edu
 *
 */


#include <assert.h>
#include "LayerPropertyROI.h"
#include "vtkRGBAColorTransferFunction.h"

using namespace std;

LayerPropertyROI::LayerPropertyROI ( QObject* parent) : LayerProperty( parent )
{
  mOpacity = 0.7;
  mRGB[0] = 1;
  mRGB[1] = 1;
  mRGB[2] = 0;

  mLUTTable = vtkSmartPointer<vtkRGBAColorTransferFunction>::New();
// mLUTTable->ClampingOff();
  this->SetColorMapChanged();

  connect( this, SIGNAL(ColorMapChanged()), this, SIGNAL(PropertyChanged()) );
  connect( this, SIGNAL(OpacityChanged(double)), this, SIGNAL(PropertyChanged()) );
}

LayerPropertyROI::~LayerPropertyROI ()
{}


vtkRGBAColorTransferFunction* LayerPropertyROI::GetLookupTable () const
{
  return mLUTTable;
}

void LayerPropertyROI::SetColorMapChanged()
{
  assert( mLUTTable.GetPointer() );

  mLUTTable->RemoveAllPoints();
  mLUTTable->AddRGBAPoint( 1-0.001, 0, 0, 0, 0 );
  mLUTTable->AddRGBAPoint( 1,    mRGB[0], mRGB[1], mRGB[2], 1 );
  mLUTTable->AddRGBAPoint( 100,  mRGB[0], mRGB[1], mRGB[2], 1 );

  mLUTTable->Build();

  // Notify the layers that use the color map stuff.
  emit ColorMapChanged();
}

void LayerPropertyROI::SetColor ( double r, double g, double b )
{
  mRGB[0] = r;
  mRGB[1] = g;
  mRGB[2] = b;
  this->SetColorMapChanged();
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

