/**
 * @file  vtkKWScubaLayer2DMRI.cxx
 * @brief A vtkKWScubaLayer that displayes 2DMRI slices
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
#include <assert.h>
#include "vtkKWScubaLayer2DMRI.h"
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

using namespace std;

vtkStandardNewMacro( vtkKWScubaLayer2DMRI );
vtkCxxRevisionMacro( vtkKWScubaLayer2DMRI, "$Revision: 1.5 $" );

vtkKWScubaLayer2DMRI::vtkKWScubaLayer2DMRI () :
  mMRIProperties( NULL )
{
  mWorldCenter[0] = mWorldCenter[1] = mWorldCenter[2] = 0;
  mWorldSize[0] = mWorldSize[1] = mWorldSize[2] = 0;
}

vtkKWScubaLayer2DMRI::~vtkKWScubaLayer2DMRI () {
}

void
vtkKWScubaLayer2DMRI::SetMRIProperties ( ScubaCollectionPropertiesMRI const* iProperties ) {
  mMRIProperties = iProperties;
}

void
vtkKWScubaLayer2DMRI::Create () {

}

void
vtkKWScubaLayer2DMRI::LoadDataFromProperties () {

  assert( mMRIProperties );

  //
  // Source object reads the volume and outputs structured points. If
  // we don't have one, we don't have anything to load.
  //
  vtkFSVolumeSource* source = mMRIProperties->GetSource();
  if( NULL == source )
    return;

  // Get some values from the MRI.
  mWorldCenter[0] = source->GetRASCenterX();
  mWorldCenter[1] = source->GetRASCenterY();
  mWorldCenter[2] = source->GetRASCenterZ();

  float RASBounds[6];
  source->GetRASBounds( RASBounds );

  mWorldSize[0] = RASBounds[1] - RASBounds[0];
  mWorldSize[1] = RASBounds[3] - RASBounds[2];
  mWorldSize[2] = RASBounds[5] - RASBounds[4];

  //
  // This transforms the voxel space source into RAS space.
  //
  vtkSmartPointer<vtkImageReslice> volumeToRAS = 
    vtkSmartPointer<vtkImageReslice>::New();
  volumeToRAS->SetInputConnection( source->GetOutputPort() );
  volumeToRAS->SetOutputDimensionality( 3 );

  // This rotates the volume to the proper orientation. From
  // ImageReslice: "applying a transform to the resampling grid (which
  // lies in the output coordinate system) is equivalent to applying
  // the inverse of that transform to the input volume."
  double* rtv = source->GetRASToVoxelMatrix();

  // 0,0 0,1 0,2 0,3      rtv[0]  rtv[1]  rtv[2]  0
  // 1,0 1,1 1,2 1,3  =>  rtv[4]  rtv[5]  rtv[6]  0
  // 2,0 2,1 2,2 2,3      rtv[8]  rtv[9]  rtv[10] 0
  // 3,0 3,1 3,2 3,3        0       0       0     1
  vtkSmartPointer<vtkMatrix4x4> matrix = 
    vtkSmartPointer<vtkMatrix4x4>::New();
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

  vtkSmartPointer<vtkTransform> transform = 
    vtkSmartPointer<vtkTransform>::New();
  transform->SetMatrix( matrix );

  volumeToRAS->SetResliceTransform( transform );
  volumeToRAS->BorderOff();

  // This sets our output extent.
  volumeToRAS->SetOutputExtent( (int)RASBounds[0], (int)RASBounds[1],
                                (int)RASBounds[2], (int)RASBounds[3],
                                (int)RASBounds[4], (int)RASBounds[5] );


  //
  // The reslice object just takes a slice out of the volume.
  //
  mReslice = vtkSmartPointer<vtkImageReslice>::New();
  mReslice->SetInputConnection( volumeToRAS->GetOutputPort() );
  mReslice->BorderOff();

  // This sets us to extract slices.
  mReslice->SetOutputDimensionality( 2 );

  // This will change depending what orienation we're in.
  mReslice->SetResliceAxesDirectionCosines( 1, 0, 0,
					    0, 1, 0,
					    0, 0, 1 );

  // This will change to select a different slice.
  mReslice->SetResliceAxesOrigin( 0, 0, 0 );

  //
  // Flip over the x axis (left/right). This get us into neurological
  // view.
  //
  vtkSmartPointer<vtkImageFlip> imageFlip = 
    vtkSmartPointer<vtkImageFlip>::New();
  imageFlip->SetInputConnection( mReslice->GetOutputPort() );
  imageFlip->SetFilteredAxis( 0 ); // x axis


  //
  // Image to colors using color table.
  //
  mColorMap = vtkSmartPointer<vtkImageMapToColors>::New();
  mColorMap->SetInputConnection( imageFlip->GetOutputPort() );
  mColorMap->SetOutputFormatToRGBA();
  mColorMap->PassAlphaToOutputOn();
  mColorMap->SetLookupTable( mMRIProperties->GetGrayScaleTable() );

  //
  // Colors to texture.
  //
  mTexture = vtkSmartPointer<vtkTexture>::New();
  mTexture->SetInputConnection( mColorMap->GetOutputPort() );
  mTexture->RepeatOff();
  mTexture->InterpolateOff();

  //
  // Plane mesh object.
  //
  vtkSmartPointer<vtkPlaneSource> plane = 
    vtkSmartPointer<vtkPlaneSource>::New();

  //
  // Plane mapper transform.
  //
  mPlaneTransform = vtkSmartPointer<vtkTransform>::New();

  //
  // Poly data from plane and plane transform.
  //
  vtkSmartPointer<vtkTransformPolyDataFilter> planePDF = 
    vtkSmartPointer<vtkTransformPolyDataFilter>::New();
  planePDF->SetInput( plane->GetOutput() );
  planePDF->SetTransform( mPlaneTransform );

  //
  // Mapper for plane.
  //
  mPlaneMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  mPlaneMapper->ImmediateModeRenderingOn();
  mPlaneMapper->SetInputConnection( planePDF->GetOutputPort() );


  //
  // Prop in scene with plane mesh and texture.
  //
  mPlaneActor = vtkSmartPointer<vtkActor>::New();
  mPlaneActor->SetMapper( mPlaneMapper );
  mPlaneActor->SetTexture( mTexture );

  // Add it to our list to render. Link it to us in the map.
  this->AddProp( mPlaneActor ); 

  // Set ourselves up.
  this->UpdateOpacity();
  this->UpdateColorMap();
  this->UpdateResliceInterpolation();
  this->UpdateTextureSmoothing();
  this->Update2DInfo();
}

void
vtkKWScubaLayer2DMRI::AddControls ( vtkKWWidget* iPanel ) {
}

void
vtkKWScubaLayer2DMRI::RemoveControls () {
}

void
vtkKWScubaLayer2DMRI::DoListenToMessage ( string const isMessage,
					  void* const iData ) {

  if( isMessage == "OpacityChanged" ) {
    this->UpdateOpacity();

  } else if( isMessage == "ColorMapChanged" ) {
    this->UpdateColorMap();

  } else if( isMessage == "ResliceInterpolationChanged" ) {
    this->UpdateResliceInterpolation();

  } else if( isMessage == "TextureSmoothingChanged" ) {
    this->UpdateTextureSmoothing();

  } else if( isMessage == "Layer2DInfoChanged" ) {
    this->Update2DInfo ();
    
  }
}

void
vtkKWScubaLayer2DMRI::GetRASBounds ( float ioBounds[6] ) const {

  assert( mMRIProperties );

  if ( mMRIProperties->GetSource() )
    mMRIProperties->GetSource()->GetRASBounds( ioBounds );
  else {
    for ( int nBound = 0; nBound < 6; nBound++ )
      ioBounds[nBound] = 0;
  }
}

void
vtkKWScubaLayer2DMRI::Get2DRASZIncrementHint ( float ioHint[3]) const {

  if ( mMRIProperties ) {
    ioHint[0] = mMRIProperties->GetSource()->GetPixelSizeX();
    ioHint[1] = mMRIProperties->GetSource()->GetPixelSizeY();
    ioHint[2] = mMRIProperties->GetSource()->GetPixelSizeZ();
  } else {
    ioHint[0] = 0;
    ioHint[1] = 0;
    ioHint[2] = 0;
  }
}

void
vtkKWScubaLayer2DMRI::GetInfoItems ( float iRAS[3],
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
vtkKWScubaLayer2DMRI::GetSource () const {
  return mMRIProperties->GetSource();
}

void
vtkKWScubaLayer2DMRI::UpdateOpacity () {

  assert( mMRIProperties );

  if ( mPlaneActor.GetPointer() )
    if ( mPlaneActor->GetProperty() ) {
      mPlaneActor->GetProperty()->SetOpacity( mProperties->GetOpacity() );
      this->PipelineChanged();
    }
}

void
vtkKWScubaLayer2DMRI::UpdateColorMap () {

  assert( mMRIProperties );

  switch ( mMRIProperties->GetColorMap() ) {
  case ScubaCollectionPropertiesMRI::NoColorMap:
    mColorMap->SetLookupTable( NULL );
    this->PipelineChanged();
    break;
    
  case ScubaCollectionPropertiesMRI::GrayScale:
    mColorMap->SetLookupTable( mMRIProperties->GetGrayScaleTable() );
    this->PipelineChanged();
    break;
    
  case ScubaCollectionPropertiesMRI::HeatScale:
    mColorMap->SetLookupTable( mMRIProperties->GetHeatScaleTable() );
    this->PipelineChanged();
    break;
    
  case ScubaCollectionPropertiesMRI::LUT:
    mColorMap->SetLookupTable( mMRIProperties->GetLUTTable() );
    this->PipelineChanged();
    break;

  default:
    break;
  }
}

void
vtkKWScubaLayer2DMRI::UpdateResliceInterpolation () {

  assert( mMRIProperties );

  if( mReslice.GetPointer() ) {
    mReslice->
      SetInterpolationMode( mMRIProperties->GetResliceInterpolation() );
    this->PipelineChanged();
  }
}

void
vtkKWScubaLayer2DMRI::UpdateTextureSmoothing () {

  assert( mMRIProperties );

  if( mTexture.GetPointer() ) {
    mTexture->SetInterpolate( mMRIProperties->GetTextureSmoothing() );
    this->PipelineChanged();
  }
}

void
vtkKWScubaLayer2DMRI::Update2DInfo () {
  
  assert( mProperties );
  assert( mViewProperties );

  float rasZ = mViewProperties->Get2DRASZ();

  switch ( mViewProperties->Get2DInPlane() ) {
  case 0:
    mPlaneTransform->Identity();
    mPlaneTransform->Translate( 0, mWorldCenter[1], mWorldCenter[2] );
    //    mPlaneTransform->Translate( mWorldCenter[1], mWorldCenter[2], 0 );
    mPlaneTransform->RotateX( 90 );
    mPlaneTransform->RotateY( 90 );
    mPlaneTransform->Scale( mWorldSize[1], mWorldSize[2], 1 );
    
    // Putting negatives in the reslice axes cosines will flip the
    // image on that axis.
    mReslice->SetResliceAxesDirectionCosines( 0, -1, 0,
						0, 0, 1,
					      1, 0, 0 );
    mReslice->SetResliceAxesOrigin( (int)(rasZ - mWorldCenter[0]), 0, 0 );
    break;
  case 1:
    mPlaneTransform->Identity();
    mPlaneTransform->Translate( mWorldCenter[0], 0, mWorldCenter[2] );
    mPlaneTransform->RotateX( 90 );
    mPlaneTransform->RotateY( 180 );
    mPlaneTransform->Scale( mWorldSize[0], mWorldSize[2], 1 );
    
    // Putting negatives in the reslice axes cosines will flip the
      // image on that axis.
    mReslice->SetResliceAxesDirectionCosines( 1, 0, 0,
					      0, 0, 1,
					      0, 1, 0 );
    mReslice->SetResliceAxesOrigin( 0, (int)(rasZ - mWorldCenter[1]), 0 );
    break;
  case 2:
    mPlaneTransform->Identity();
    mPlaneTransform->Translate( mWorldCenter[0], mWorldCenter[1], 0 );
    mPlaneTransform->RotateY( 180 );
    mPlaneTransform->Scale( mWorldSize[0], mWorldSize[1], 1 );
    
    mReslice->SetResliceAxesDirectionCosines( 1, 0, 0,
					      0, 1, 0,
					      0, 0, 1 );
    mReslice->SetResliceAxesOrigin( 0, 0, (int)(rasZ - mWorldCenter[2]) );
    break;
  }
}

