/**
 * @file  vtkKWScubaLayer2DMRIS.cxx
 * @brief Displays 2D slices of surfaces
 *
 * Displays surfaces where they intersect with a 2D slice.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: kteich $
 *    $Date: 2007/04/06 22:23:04 $
 *    $Revision: 1.1 $
 *
 * Copyright (C) 2002-2007,
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


#include <string>
#include <stdexcept>
#include "vtkKWScubaLayer2DMRIS.h"
#include "ScubaCollectionProperties.h"
#include "ScubaCollectionPropertiesMRIS.h"
#include "vtkFSSurfaceSource.h"
#include "vtkObjectFactory.h"
#include "vtkPolyDataMapper.h"
#include "vtkActor.h"
#include "vtkDecimatePro.h"
#include "vtkPlane.h"
#include "vtkCutter.h"
#include "vtkProperty.h"
#include "vtkTransformPolyDataFilter.h"
#include "vtkMatrix4x4.h"
#include "vtkTransform.h"
#include "vtkKWRadioButtonSet.h"
#include "vtkKWRadioButton.h"
#include "vtkKWLabel.h"
#include "ScubaInfoItem.h"

using namespace std;

vtkStandardNewMacro( vtkKWScubaLayer2DMRIS );
vtkCxxRevisionMacro( vtkKWScubaLayer2DMRIS, "$Revision: 1.1 $" );

vtkKWScubaLayer2DMRIS::vtkKWScubaLayer2DMRIS () :
  mMRISProperties( NULL ),
  mSlicePlane( NULL ),
  mNormalMapper( NULL ),
  mFastMapper( NULL ),
  mActor( NULL ) {
  mRASCenter[0] = mRASCenter[1] = mRASCenter[2] = 0;
}

vtkKWScubaLayer2DMRIS::~vtkKWScubaLayer2DMRIS () {

  if ( mNormalMapper )
    mNormalMapper->Delete();
  if ( mFastMapper )
    mFastMapper->Delete();
  if ( mSlicePlane )
    mSlicePlane->Delete();
  if ( mActor )
    mActor->Delete();
}

void
vtkKWScubaLayer2DMRIS::SetMRISProperties ( ScubaCollectionPropertiesMRIS* const iProperties ) {
  mMRISProperties = iProperties;
}


void
vtkKWScubaLayer2DMRIS::Create () {

  // Bail if we don't have our source yet.
  if( NULL == mMRISProperties )
    throw runtime_error( "vtkKWScubaLayer2DMRIS::Create: No source" );
  
  //
  // Source object reads the surface and outputs poly data.
  //
  vtkFSSurfaceSource* source = mMRISProperties->GetSource();

  // Get some info from the MRIS.
  mRASCenter[0] = source->GetRASCenterX();
  mRASCenter[1] = source->GetRASCenterY();
  mRASCenter[2] = source->GetRASCenterZ();

  //
  // Cutting plane.
  //
  float rasZ = mViewProperties->Get2DRASZ();
  int inPlane = mViewProperties->Get2DInPlane();
  mSlicePlane = vtkPlane::New();
  mSlicePlane->SetOrigin( rasZ, rasZ, rasZ );
  mSlicePlane->SetNormal( (inPlane==0), (inPlane==1), (inPlane==2) );

  //
  // Cutters that takes the 3D surface and outputs 2D lines.
  //
  vtkCutter* clipper = vtkCutter::New();
  clipper->SetInputConnection( mMRISProperties->GetNormalModeOutputPort() );
  clipper->SetCutFunction( mSlicePlane );

  vtkCutter* fastClipper = vtkCutter::New();
  fastClipper->SetInputConnection( mMRISProperties->GetFastModeOutputPort() );
  fastClipper->SetCutFunction( mSlicePlane );

  //
  // Mappers for the lines.
  //
  mNormalMapper = vtkPolyDataMapper::New();
  mNormalMapper->SetInput( clipper->GetOutput() );

  mFastMapper = vtkPolyDataMapper::New();
  mFastMapper->SetInput( fastClipper->GetOutput() );

  //
  // Actors in the scene, drawing the mapped lines.
  //
  mActor = vtkActor::New();
  mActor->SetMapper( mNormalMapper );
  mActor->SetBackfaceProperty( mActor->MakeProperty() );
  mActor->GetBackfaceProperty()->BackfaceCullingOff();
  mActor->SetProperty( mActor->MakeProperty() );
  mActor->GetProperty()->SetColor( 1, 0, 0 );
  mActor->GetProperty()->SetEdgeColor( 1, 0, 0 );
  mActor->GetProperty()->SetInterpolationToFlat();

  this->AddProp( mActor );

  
  // Update the fastClipper now so we do the decimation up front.
  mFastMapper->Update();
}

void
vtkKWScubaLayer2DMRIS::AddControls ( vtkKWWidget* iPanel ) {

}

void
vtkKWScubaLayer2DMRIS::RemoveControls () {

}

void
vtkKWScubaLayer2DMRIS::DoListenToMessage ( string isMessage, void* iData ) {

  if( isMessage == "OpacityChanged" ) {

    if ( mActor )
      if ( mActor->GetProperty() ) {
	mActor->GetProperty()->SetOpacity( mProperties->GetOpacity() );
	this->PipelineChanged();
      }

  } else if( isMessage == "FastModeChanged" ) {
    
    if ( mViewProperties->GetFastMode() )
      mActor->SetMapper( mFastMapper );
    else
      mActor->SetMapper( mNormalMapper );
    
    this->PipelineChanged();

  } else if( isMessage == "Layer2DInfoChanged" ) {

    if ( mSlicePlane ) {

      float rasZ = mViewProperties->Get2DRASZ();
      int inPlane = mViewProperties->Get2DInPlane();
    
      // Update the location of the slice plane. We also transform the
      // RAS by the world center here.
      mSlicePlane->SetOrigin( (inPlane==0) ? rasZ - mRASCenter[0] : 0,
			      (inPlane==1) ? rasZ - mRASCenter[1] : 0,
			      (inPlane==2) ? rasZ - mRASCenter[2] : 0);
      mSlicePlane->SetNormal( (inPlane==0), (inPlane==1), (inPlane==2) );
      
      // When we slice the polydata, the physical location of the lines
      // in 3D space is that of the RAS plane on which we sliced, so we
      // need to transform the resulting data onto the 0 plane. However,
      // for some weird reason, the surface will disappear under the
      // texture planes in some slices, so move it just a little in
      // front as well. Also adjust by the world center.
      mActor->SetPosition( (inPlane==0) ? -(rasZ - mRASCenter[0]) + 0.1 :
			   mRASCenter[0],
			   (inPlane==1) ? -(rasZ - mRASCenter[1]) + 0.1 :
			   mRASCenter[1],
			   (inPlane==2) ? -(rasZ - mRASCenter[2]) - 0.1 :
			   mRASCenter[2] );
      
      this->PipelineChanged();
    }
  }
}

void
vtkKWScubaLayer2DMRIS::GetRASBounds ( float ioBounds[6] ) const {

  if ( mMRISProperties )
    mMRISProperties->GetSource()->GetRASBounds( ioBounds );
  else {
    for ( int nBound = 0; nBound < 6; nBound++ )
      ioBounds[nBound] = 0;
  }
}

void
vtkKWScubaLayer2DMRIS::Get2DRASZIncrementHint ( float ioHint[3]) const {

  ioHint[0] = 1.0;
  ioHint[1] = 1.0;
  ioHint[2] = 1.0;
}

void
vtkKWScubaLayer2DMRIS::GetInfoItems ( float iRAS[3],
                                      list<ScubaInfoItem>& ilInfo ) const {

  ScubaInfoItem info;
  char sLabel[1024], sValue[1024];

  try {

    float distance;
    int nVertex =
      mMRISProperties->GetSource()->FindVertexAtSurfaceRAS( iRAS, &distance );

    sprintf( sLabel, "%s vno", mProperties->GetLabel() );
    sprintf( sValue, "%d (%.2f)", nVertex, distance );
    info.Clear();
    info.SetLabel( sLabel );
    info.SetValue( sValue );
    info.SetShortenHint( false );
  } catch (...) {

    // No vertex found, don't worry.
    sprintf( sLabel, "%s vno", mProperties->GetLabel() );
    strncpy( sValue, "NONE", sizeof(sValue) );
    info.Clear();
    info.SetLabel( sLabel );
    info.SetValue( sValue );
    info.SetShortenHint( false );
  }

  // Return it.
  ilInfo.push_back( info );
}

