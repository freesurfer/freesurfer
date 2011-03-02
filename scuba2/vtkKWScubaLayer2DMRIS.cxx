/**
 * @file  vtkKWScubaLayer2DMRIS.cxx
 * @brief Displays 2D slices of surfaces
 *
 * Displays surfaces where they intersect with a 2D slice.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:40 $
 *    $Revision: 1.6 $
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


#include <string>
#include <stdexcept>
#include <assert.h>
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
#include "vtkPointData.h"
#include "vtkTransform.h"
#include "vtkKWRadioButtonSet.h"
#include "vtkKWRadioButton.h"
#include "vtkKWLabel.h"
#include "ScubaInfoItem.h"

using namespace std;

vtkStandardNewMacro( vtkKWScubaLayer2DMRIS );
vtkCxxRevisionMacro( vtkKWScubaLayer2DMRIS, "$Revision: 1.6 $" );

vtkKWScubaLayer2DMRIS::vtkKWScubaLayer2DMRIS () :
  mMRISProperties( NULL ) {
}

vtkKWScubaLayer2DMRIS::~vtkKWScubaLayer2DMRIS () {
}

void
vtkKWScubaLayer2DMRIS::SetMRISProperties ( ScubaCollectionPropertiesMRIS const* iProperties ) {
  mMRISProperties = iProperties;
}


void
vtkKWScubaLayer2DMRIS::Create () {

}

void
vtkKWScubaLayer2DMRIS::LoadDataFromProperties () {

  assert( mMRISProperties );
  assert( mViewProperties );

  //
  // Source object reads the surface and outputs poly data. If we
  // don't have one, we don't have anything to load.
  //
  vtkFSSurfaceSource* source = mMRISProperties->GetSource();
  if( NULL == source )
    return;

  //
  // Cutting plane.
  //
  float rasZ = mViewProperties->Get2DRASZ();
  int inPlane = mViewProperties->Get2DInPlane();
  mSlicePlane = vtkSmartPointer<vtkPlane>::New();
  mSlicePlane->SetOrigin( rasZ, rasZ, rasZ );
  mSlicePlane->SetNormal( (inPlane==0), (inPlane==1), (inPlane==2) );

  //
  // Cutters that takes the 3D surface and outputs 2D lines.
  //
  vtkSmartPointer<vtkCutter> clipper = 
    vtkSmartPointer<vtkCutter>::New();
  clipper->SetInputConnection( mMRISProperties->GetNormalModeOutputPort() );
  clipper->SetCutFunction( mSlicePlane );

  vtkSmartPointer<vtkCutter> fastClipper = 
    vtkSmartPointer<vtkCutter>::New();
  fastClipper->SetInputConnection( mMRISProperties->GetFastModeOutputPort() );
  fastClipper->SetCutFunction( mSlicePlane );

  //
  // Mappers for the lines.
  //
  mNormalMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  mNormalMapper->SetInput( clipper->GetOutput() );

  mFastMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  mFastMapper->SetInput( fastClipper->GetOutput() );

  //
  // Actors in the scene, drawing the mapped lines.
  //
  mActor = vtkSmartPointer<vtkActor>::New();
  mActor->SetMapper( mNormalMapper );
  vtkSmartPointer<vtkProperty> property;
  property.TakeReference( mActor->MakeProperty() );
  mActor->SetBackfaceProperty( property );
  property->BackfaceCullingOff();
  property.TakeReference( mActor->MakeProperty() );
  mActor->SetProperty( property );
  property->SetColor( 1, 0, 0 );
  property->SetEdgeColor( 1, 0, 0 );
  property->SetInterpolationToFlat();

  this->AddProp( mActor );

  
  // Update the fastClipper now so we do the decimation up front.
  mFastMapper->Update();

  // If we have scalars, display thouse.
  if( mMRISProperties->GetScalarsValues() ) {

    vtkFloatArray* scalars = mMRISProperties->GetScalarsValues();
    mNormalMapper->SetColorModeToMapScalars();
    mFastMapper->SetColorModeToMapScalars();

    source->GetOutput()->GetPointData()->
      SetScalars( reinterpret_cast<vtkDataArray*>(scalars) );

    if( mMRISProperties->GetScalarsColors() ) {
      
      vtkScalarsToColors* colors = mMRISProperties->GetScalarsColors();
      mNormalMapper->UseLookupTableScalarRangeOn();
      mNormalMapper->SetLookupTable( colors );
      mFastMapper->UseLookupTableScalarRangeOn();
      mFastMapper->SetLookupTable( colors );
    }

  }
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
    
      // Update the location of the slice place to extract the proper
      // slice.
      mSlicePlane->SetOrigin( (inPlane==0) ? rasZ : 0,
			      (inPlane==1) ? rasZ : 0,
			      (inPlane==2) ? rasZ : 0);

      // Set the plane's orientation.
      mSlicePlane->SetNormal( (inPlane==0), (inPlane==1), (inPlane==2) );

      // When we slice the polydata, the physical location of the lines
      // in 3D space is that of the RAS plane on which we sliced, so we
      // need to transform the resulting data onto the 0 plane. However,
      // for some weird reason, the surface will disappear under the
      // texture planes in some slices, so move it just a little in
      // front as well.
      mActor->SetPosition( (inPlane==0) ? -rasZ + 0.1 : 0,
			   (inPlane==1) ? -rasZ + 0.1 : 0,
			   (inPlane==2) ? -rasZ - 0.1 : 0);
 
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
      mMRISProperties->GetSource()->FindVertexAtRAS( iRAS, &distance );

    sprintf( sLabel, "%s vno", mProperties->GetLabel() );
    double foundPointRAS[3];
    mMRISProperties->GetSource()->GetSurfaceRASAtVertex(nVertex,foundPointRAS);
    sprintf( sValue, "%d (%.2f) %.2f %.2f %.2f", nVertex, distance,
	     (float)foundPointRAS[0], (float)foundPointRAS[1], (float)foundPointRAS[2]);
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

