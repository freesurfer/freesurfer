/**
 * @file  vtkKWScubaLayer3DMRI.cxx
 * @brief A vtkKWScubaLayer that displayes 3DMRI slices
 *
 * Displayes MRI volumes (from a vtkFSVolumeSource) in a slice.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:40 $
 *    $Revision: 1.5 $
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
#include "vtkKWScubaLayer3DMRI.h"
#include "ScubaCollectionProperties.h"
#include "ScubaCollectionPropertiesMRI.h"
#include "vtkObjectFactory.h"
#include "vtkFSVolumeSource.h"
#include "vtkImageReslice.h"
#include "vtkImageMapToColors.h"
#include "vtkTransform.h"
#include "vtkTexture.h"
#include "vtkPolyDataMapper.h"
#include "vtkActor.h"
#include "vtkMatrix4x4.h"
#include "vtkImageFlip.h"
#include "vtkPlaneSource.h"
#include "vtkTransformPolyDataFilter.h"
#include "vtkFreesurferLookupTable.h"
#include "vtkRGBATransferFunction.h"
#include "vtkProperty.h"
#include "vtkKWScaleWithEntry.h"

using namespace std;

vtkStandardNewMacro( vtkKWScubaLayer3DMRI );
vtkCxxRevisionMacro( vtkKWScubaLayer3DMRI, "$Revision: 1.5 $" );

vtkKWScubaLayer3DMRI::vtkKWScubaLayer3DMRI () :
  mMRIProperties( NULL )
{
  for( int n = 0; n < 3; n++ ) {
    mReslice[n] = NULL;
    mColorMap[n] = NULL;
    mPlaneTransform[n] = NULL;
    mTexture[n] = NULL;
    mPlaneMapper[n] = NULL;
    mPlaneActor[n] = NULL;
    mWorldCenter[n] = 0;
    mWorldSize[n] = 0;
  }
}

vtkKWScubaLayer3DMRI::~vtkKWScubaLayer3DMRI () {
}

void
vtkKWScubaLayer3DMRI::SetMRIProperties ( ScubaCollectionPropertiesMRI const* iProperties ) {
  mMRIProperties = iProperties;
}

void
vtkKWScubaLayer3DMRI::Create () {

  // Bail if we don't have our source and tables yet.
  if( NULL == mMRIProperties )
    throw runtime_error( "vtkKWScubaLayer3DMRI::Create: No source" );
  
  //
  // Source object reads the volume and outputs structured points.
  //
  vtkFSVolumeSource* source = mMRIProperties->GetSource();

  // Get some values from the MRI.
  mWorldCenter[0] = source->GetRASCenterX();
  mWorldCenter[1] = source->GetRASCenterY();
  mWorldCenter[2] = source->GetRASCenterZ();

  float RASBounds[6];
  this->GetRASBounds( RASBounds );

  mWorldSize[0] = RASBounds[1] - RASBounds[0];
  mWorldSize[1] = RASBounds[3] - RASBounds[2];
  mWorldSize[2] = RASBounds[5] - RASBounds[4];

  //
  // This transforms the voxel space source into RAS space.
  //
  vtkImageReslice* volumeToRAS = vtkImageReslice::New();
  volumeToRAS->SetInputConnection( source->GetOutputPort() );
  volumeToRAS->SetOutputDimensionality( 3 );

  // This rotates the volume to the proper orientation. From
  // ImageReslice: "applying a transform to the resampling grid (which
  // lies in the output coordinate system) is equivalent to applying
  // the inverse of that transform to the input volume."
//  double rtv[ 16 ];
//  source->GetUnscaledRASToVoxelMatrix( rtv );

  double *rtv = source->GetRASToVoxelMatrix();

  vtkMatrix4x4* matrix = vtkMatrix4x4::New();
  matrix->SetElement( 0, 0, rtv[0] );
  matrix->SetElement( 0, 1, rtv[1] );
  matrix->SetElement( 0, 2, rtv[2] );
  matrix->SetElement( 0, 3, 0 );
  matrix->SetElement( 1, 0, rtv[4] );
  matrix->SetElement( 1, 1, rtv[5] );
  matrix->SetElement( 1, 2, rtv[6] );
  matrix->SetElement( 1, 3, 0 );
  matrix->SetElement( 2, 0, rtv[8] );
  matrix->SetElement( 2, 1, rtv[9] );
  matrix->SetElement( 2, 2, rtv[10] );
  matrix->SetElement( 2, 3, 0 );
  matrix->SetElement( 3, 0, 0 );
  matrix->SetElement( 3, 1, 0 );
  matrix->SetElement( 3, 2, 0 );
  matrix->SetElement( 3, 3, 1 );

  vtkTransform* transform = vtkTransform::New();
  transform->SetMatrix( matrix );
  matrix->Delete();
  
  volumeToRAS->SetResliceTransform( transform );
  volumeToRAS->BorderOff();
  transform->Delete();

  // This sets our output extent.
  volumeToRAS->SetOutputExtent( (int)RASBounds[0], (int)RASBounds[1],
                                (int)RASBounds[2], (int)RASBounds[3],
                                (int)RASBounds[4], (int)RASBounds[5] );


  // This is like the code in 2DMRI, but we do it three times, since
  // we're creating three planes.
  for( int n = 0; n < 3; n++ ) {

    //
    // The reslice object just takes a slice out of the volume.
    //
    if ( !mReslice[n] )
      mReslice[n] = vtkImageReslice::New();
    mReslice[n]->SetInputConnection( volumeToRAS->GetOutputPort() );
    mReslice[n]->BorderOff();
    
    // This sets us to extract slices.
    mReslice[n]->SetOutputDimensionality( 2 );
    
    // This will change depending what orienation we're in.
    mReslice[n]->SetResliceAxesDirectionCosines( 1, 0, 0,
						 0, 1, 0,
						 0, 0, 1 );
    
    // This will change to select a different slice.
    mReslice[n]->SetResliceAxesOrigin( 0, 0, 0 );
    
    //
    // Flip over the x axis (left/right). This get us into neurological
    // view.
    //
    vtkImageFlip* imageFlip = vtkImageFlip::New();
    imageFlip->SetInputConnection( mReslice[n]->GetOutputPort() );
    imageFlip->SetFilteredAxis( 0 ); // x axis
    
    //
    // Image to colors using color table.
    //
    mColorMap[n] = vtkImageMapToColors::New();
    mColorMap[n]->SetInputConnection( imageFlip->GetOutputPort() );
    mColorMap[n]->SetOutputFormatToRGBA();
    mColorMap[n]->PassAlphaToOutputOn();
    mColorMap[n]->SetLookupTable( mMRIProperties->GetGrayScaleTable() );
    imageFlip->Delete();
    
    //
    // Colors to texture.
    //
    if ( !mTexture[n] )
      mTexture[n] = vtkTexture::New();
    mTexture[n]->SetInputConnection( mColorMap[n]->GetOutputPort() );
    mTexture[n]->RepeatOff();
    mTexture[n]->InterpolateOff();
    
    //
    // Plane mesh object.
    //
    vtkPlaneSource* plane = vtkPlaneSource::New();
    
    //
    // Plane mapper transform.
    //
    if ( !mPlaneTransform[n] )
      mPlaneTransform[n] = vtkTransform::New();
    
    //
    // Poly data from plane and plane transform.
    //
    vtkTransformPolyDataFilter* planePDF = vtkTransformPolyDataFilter::New();
    planePDF->SetInput( plane->GetOutput() );
    planePDF->SetTransform( mPlaneTransform[n] );
    plane->Delete();
    
    //
    // Mapper for plane.
    //
    mPlaneMapper[n] = vtkPolyDataMapper::New();
    mPlaneMapper[n]->ImmediateModeRenderingOn();
    mPlaneMapper[n]->SetInputConnection( planePDF->GetOutputPort() );
    planePDF->Delete();

    
    //
    // Prop in scene with plane mesh and texture.
    //
    if ( !mPlaneActor[n] )
      mPlaneActor[n] = vtkActor::New();
    mPlaneActor[n]->SetMapper( mPlaneMapper[n] );
    mPlaneActor[n]->SetTexture( mTexture[n] );

    mPlaneActor[n]->GetProperty()->SetColor( 1.0, 1.0, 1.0 );
    mPlaneActor[n]->GetProperty()->SetAmbient( 1 );
    
    // Add it to our list to render. Link it to us in the map.
    this->AddProp( mPlaneActor[n] ); 
  }

  // Delete the stuff now.
  volumeToRAS->Delete();

  // Set us up with stuff from the properties.
  this->UpdateOpacity();
  this->UpdateColorMap();
  this->UpdateResliceInterpolation();
  this->UpdateTextureSmoothing();
  this->UpdatePlanes();
}

void
vtkKWScubaLayer3DMRI::DoListenToMessage ( string const isMessage,
					  void* const iData ) {

  if( isMessage == "OpacityChanged" ) {
    this->UpdateOpacity();
    
  } else if( isMessage == "ColorMapChanged" ) {
    this->UpdateColorMap();

  } else if( isMessage == "ResliceInterpolationChanged" ) {
    this->UpdateResliceInterpolation();
    
  } else if( isMessage == "TextureSmoothingChanged" ) {
    this->UpdateTextureSmoothing();

  } else if( isMessage == "Layer3DInfoChanged" ) {
    this->UpdatePlanes();

  }
}

void
vtkKWScubaLayer3DMRI::GetRASBounds ( float ioBounds[6] ) const {

  if ( mMRIProperties )
    mMRIProperties->GetSource()->GetRASBounds( ioBounds );
  else {
    for ( int nBound = 0; nBound < 6; nBound++ )
      ioBounds[nBound] = 0;
  }
}

void
vtkKWScubaLayer3DMRI::GetInfoItems ( float iRAS[3],
                                     list<ScubaInfoItem>& ilInfo ) const {

  ScubaInfoItem info;

  int idx[3];
  vtkFSVolumeSource* source = mMRIProperties->GetSource();
  source->ConvertRASToIndex( iRAS[0], iRAS[1], iRAS[2],
			     idx[0], idx[1], idx[2] );

  // Build our current info item.
  char sLabel[1024];
  sprintf( sLabel, "%s index", mProperties->GetLabel() );
  info.Clear();
  info.SetLabel( sLabel );

  char sIdx[1024];
  snprintf( sIdx, sizeof(sIdx), "%d %d %d", idx[0], idx[1], idx[2] );
  info.SetValue( sIdx );
  info.SetShortenHint( false );

  // Return the info.
  ilInfo.push_back( info );
  
  // If we have an LUT color table, return the label associated with
  // the value here, otherwise just return the value.
  sprintf( sLabel, "%s value", mProperties->GetLabel() );
  info.Clear();
  info.SetLabel( sLabel );
  char sValue[1024];
  if ( idx[0] >= 0 && idx[0] < source->GetXDimension() &&
       idx[1] >= 0 && idx[1] < source->GetYDimension() &&
       idx[2] >= 0 && idx[2] < source->GetZDimension() ) {
    if( ScubaCollectionPropertiesMRI::LUT == mMRIProperties->GetColorMap() && 
    NULL != mMRIProperties->GetLUTCTAB() ) {
      int nEntry = (int)source->GetValueAtIndex( idx[0], idx[1], idx[2] );
      strncpy( sValue, "None", sizeof(sValue) );
      CTABcopyName( mMRIProperties->GetLUTCTAB(), 
		    nEntry, sValue, sizeof(sValue) );
    } else {
      snprintf( sValue, sizeof(sValue), "%.2f",
		source->GetValueAtIndex( idx[0], idx[1], idx[2] ) );
    }
  } else {
    strncpy( sValue, "OOB", sizeof(sValue) );
  }
  info.SetValue( sValue );
  info.SetShortenHint( false );

  // Return the info.
  ilInfo.push_back( info );
}

vtkFSVolumeSource* 
vtkKWScubaLayer3DMRI::GetSource () const {
  return mMRIProperties->GetSource();
}

void
vtkKWScubaLayer3DMRI::AddControls ( vtkKWWidget* iPanel ) {

}

void
vtkKWScubaLayer3DMRI::RemoveControls () {

}

void
vtkKWScubaLayer3DMRI::UpdateOpacity () {
  
  if( NULL == mProperties )
    return;

  for( int n = 0; n < 3; n++ )
    if ( mPlaneActor[n] )
      if ( mPlaneActor[n]->GetProperty() )
	mPlaneActor[n]->GetProperty()->
	  SetOpacity( mProperties->GetOpacity() );
  this->PipelineChanged();
}

void
vtkKWScubaLayer3DMRI::UpdateColorMap () {
  
  if( NULL == mMRIProperties )
    return;

  switch ( mMRIProperties->GetColorMap() ) {
  case ScubaCollectionPropertiesMRI::NoColorMap:
    for( int n = 0; n < 3; n++ )
      mColorMap[n]->SetLookupTable( NULL );
    this->PipelineChanged();
    break;
    
  case ScubaCollectionPropertiesMRI::GrayScale:
    for( int n = 0; n < 3; n++ )
      mColorMap[n]->SetLookupTable( mMRIProperties->GetGrayScaleTable() );
    this->PipelineChanged();
    break;
    
  case ScubaCollectionPropertiesMRI::HeatScale:
    for( int n = 0; n < 3; n++ )
      mColorMap[n]->SetLookupTable( mMRIProperties->GetHeatScaleTable() );
    this->PipelineChanged();
    break;
    
  case ScubaCollectionPropertiesMRI::LUT:
    for( int n = 0; n < 3; n++ )
      mColorMap[n]->SetLookupTable( mMRIProperties->GetLUTTable() );
    this->PipelineChanged();
    break;

  default:
    break;
  }
  
}

void
vtkKWScubaLayer3DMRI::UpdateResliceInterpolation () {  
  
  if( NULL == mMRIProperties ) 
    return;

  for( int n = 0; n < 3; n++ )
    if( mReslice[n] )
      mReslice[n]->
	SetInterpolationMode( mMRIProperties->GetResliceInterpolation() );
  this->PipelineChanged();
}

void
vtkKWScubaLayer3DMRI::UpdateTextureSmoothing () {
  
  if( NULL == mMRIProperties )
    return;
  
  for( int n = 0; n < 3; n++ )
    if( mTexture[n] )
      mTexture[n]->SetInterpolate( mMRIProperties->GetTextureSmoothing() );
  this->PipelineChanged();
}


void
vtkKWScubaLayer3DMRI::UpdatePlanes () {

  float rasX = mViewProperties->Get3DRASX();
  float rasY = mViewProperties->Get3DRASY();
  float rasZ = mViewProperties->Get3DRASZ();
  
  //
  // Now we'll set up our three planes with the right transforms.
  //
  mPlaneTransform[0]->Identity();
  mPlaneTransform[0]->Translate( (rasX - mWorldCenter[0]), 0, 0 );
  mPlaneTransform[0]->RotateX( 90 );
  mPlaneTransform[0]->RotateY( 90 );
  mPlaneTransform[0]->Scale( mWorldSize[1], mWorldSize[2], 1 );
  
  // Putting negatives in the reslice axes cosines will flip the
  // image on that axis.
  mReslice[0]->SetResliceAxesDirectionCosines( 0, -1, 0,
					       0, 0, 1,
					       1, 0, 0 );
  mReslice[0]->
    SetResliceAxesOrigin( (int)(rasX - mWorldCenter[0]), 0, 0 );
  
  
  mPlaneTransform[1]->Identity();
  mPlaneTransform[1]->Translate( 0, (rasY - mWorldCenter[1]), 0 );
  mPlaneTransform[1]->RotateX( 90 );
  mPlaneTransform[1]->RotateY( 180 );
  mPlaneTransform[1]->Scale( mWorldSize[0], mWorldSize[2], 1 );
  
  // Putting negatives in the reslice axes cosines will flip the
  // image on that axis.
  mReslice[1]->SetResliceAxesDirectionCosines( 1, 0, 0,
					       0, 0, 1,
					       0, 1, 0 );
  mReslice[1]->
    SetResliceAxesOrigin( 0, (int)(rasY - mWorldCenter[1]), 0 );
  
  
  mPlaneTransform[2]->Identity();
  mPlaneTransform[2]->Translate( 0, 0, (rasZ - mWorldCenter[2]) );
  mPlaneTransform[2]->RotateY( 180 );
  mPlaneTransform[2]->Scale( mWorldSize[0], mWorldSize[1], 1 );
  
  mReslice[2]->SetResliceAxesDirectionCosines( 1, 0, 0,
					       0, 1, 0,
					       0, 0, 1 );
  mReslice[2]->
    SetResliceAxesOrigin( 0, 0, (int)(rasZ - mWorldCenter[2]) );

  this->PipelineChanged();
}
