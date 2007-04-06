/**
 * @file  vtkKWScubaLayer3DMRIS.cxx
 * @brief Displays 3D surfaces
 *
 * Displays surfaces as a 3D model.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: kteich $
 *    $Date: 2007/04/06 22:23:05 $
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
#include "vtkKWScubaLayer3DMRIS.h"
#include "ScubaCollectionProperties.h"
#include "ScubaCollectionPropertiesMRIS.h"
#include "vtkFSSurfaceSource.h"
#include "vtkObjectFactory.h"
#include "vtkPolyDataMapper.h"
#include "vtkActor.h"
#include "vtkDecimatePro.h"
#include "vtkProperty.h"
#include "vtkTransformPolyDataFilter.h"
#include "vtkMatrix4x4.h"
#include "vtkTransform.h"
#include "vtkKWLabel.h"
#include "ScubaInfoItem.h"

using namespace std;

vtkStandardNewMacro( vtkKWScubaLayer3DMRIS );
vtkCxxRevisionMacro( vtkKWScubaLayer3DMRIS, "$Revision: 1.1 $" );

vtkKWScubaLayer3DMRIS::vtkKWScubaLayer3DMRIS () :
  mMRISProperties( NULL ),
  mNormalMapper( NULL ),
  mFastMapper( NULL ),
  mActor( NULL ) {
  mRASCenter[0] = mRASCenter[1] = mRASCenter[2] = 0;
}

vtkKWScubaLayer3DMRIS::~vtkKWScubaLayer3DMRIS () {

  if ( mNormalMapper )
    mNormalMapper->Delete();
  if ( mFastMapper )
    mFastMapper->Delete();
  if ( mActor )
    mActor->Delete();
}

void
vtkKWScubaLayer3DMRIS::SetMRISProperties ( ScubaCollectionPropertiesMRIS* const iProperties ) {
  mMRISProperties = iProperties;
}


void
vtkKWScubaLayer3DMRIS::Create () {

  // Bail if we don't have our source yet.
  if( NULL == mMRISProperties )
    throw runtime_error( "vtkKWScubaLayer3DMRIS::Create: No source" );
  
  //
  // Source object reads the surface and outputs poly data.
  //
  vtkFSSurfaceSource* source = mMRISProperties->GetSource();

  // Get some info from the MRIS.
  mRASCenter[0] = source->GetRASCenterX();
  mRASCenter[1] = source->GetRASCenterY();
  mRASCenter[2] = source->GetRASCenterZ();

  //
  // Mappers for the lines.
  //
  mNormalMapper = vtkPolyDataMapper::New();
  mNormalMapper->SetInput( mMRISProperties->GetNormalModeOutput() );

  mFastMapper = vtkPolyDataMapper::New();
  mFastMapper->SetInput( mMRISProperties->GetFastModeOutput() );

  //
  // Actors in the scene, drawing the mapped lines.
  //
  mActor = vtkActor::New();
  mActor->SetMapper( mNormalMapper );

  this->AddProp( mActor );

  
  // Update the fastClipper now so we do the decimation up front.
  mFastMapper->Update();
}

void
vtkKWScubaLayer3DMRIS::AddControls ( vtkKWWidget* iPanel ) {

}

void
vtkKWScubaLayer3DMRIS::RemoveControls () {

}

void
vtkKWScubaLayer3DMRIS::DoListenToMessage ( string const isMessage,
					   void* const iData ) {

  if( isMessage == "OpacityChanged" ) {

    if ( mActor )
      if ( mActor->GetProperty() ) {
	mActor->GetProperty()->SetOpacity( mProperties->GetOpacity() );
	this->PipelineChanged();
      }

  } else if( isMessage == "FastModeChanged" ) {
    
//     if ( mProperties->GetFastMode() )
//       mActor->SetMapper( mFastMapper );
//     else
//       mActor->SetMapper( mNormalMapper );
    
//     this->PipelineChanged();

  }
}

void
vtkKWScubaLayer3DMRIS::GetRASBounds ( float ioBounds[6] ) const {

  if ( mMRISProperties )
    mMRISProperties->GetSource()->GetRASBounds( ioBounds );
  else {
    for ( int nBound = 0; nBound < 6; nBound++ )
      ioBounds[nBound] = 0;
  }
}

void
vtkKWScubaLayer3DMRIS::GetInfoItems ( float iRAS[3],
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

