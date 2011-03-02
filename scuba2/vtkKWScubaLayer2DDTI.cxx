/**
 * @file  vtkKWScubaLayer2DDTI.cxx
 * @brief Displays 2D slices of surfaces
 *
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:40 $
 *    $Revision: 1.11 $
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

#include "vtkKWScubaLayer2DDTI.h"

#include "ScubaCollectionProperties.h"
#include "ScubaCollectionPropertiesDTI.h"
#include "ScubaInfoItem.h"
#include "vtkActor.h"
#include "vtkFDTensorGlyph.h"
#include "vtkFSVolumeSource.h"
#include "vtkImageAppendComponents.h"
#include "vtkImageClip.h"
#include "vtkImageReslice.h"
#include "vtkLODActor.h"
#include "vtkMatrix4x4.h"
#include "vtkObjectFactory.h"
#include "vtkOutlineFilter.h"
#include "vtkPolyDataMapper.h"
#include "vtkProperty.h"
#include "vtkTransform.h"
#include "vtkTransformPolyDataFilter.h"
#include "vtkImageFlip.h"

using namespace std;

vtkStandardNewMacro( vtkKWScubaLayer2DDTI );
vtkCxxRevisionMacro( vtkKWScubaLayer2DDTI, "$Revision: 1.11 $" );

vtkKWScubaLayer2DDTI::vtkKWScubaLayer2DDTI () :
  mDTIProperties( NULL ),
  mVolumeToRAS( NULL ),
  mVolumeToRASSlice( NULL ),
  mGlyph( NULL ),
  mMapper( NULL ),
  mActor( NULL ),
  mEdgeActor( NULL ),
  mPlaneTransform( NULL )
{
  for( int n = 0; n < 3; n++ ) {
    mWorldCenter[n] = 0;
    mWorldSize[n] = 0;
  }
}

vtkKWScubaLayer2DDTI::~vtkKWScubaLayer2DDTI () {

  if( mVolumeToRAS )
    mVolumeToRAS->Delete();
    
  if( mVolumeToRASSlice )
    mVolumeToRASSlice->Delete();  

  if( mGlyph )
    mGlyph->Delete();

  if( mMapper )
    mMapper->Delete();

  if ( mActor )
    mActor->Delete();

  if ( mEdgeActor )
    mEdgeActor->Delete();
    
  if( mPlaneTransform )
    mPlaneTransform->Delete();

}

void
vtkKWScubaLayer2DDTI::SetDTIProperties ( ScubaCollectionPropertiesDTI const* iProperties ) {
  mDTIProperties = iProperties;
}


void
vtkKWScubaLayer2DDTI::Create () {

  // Bail if we don't have our source yet.
  if( NULL == mDTIProperties )
    throw runtime_error( "vtkKWScubaLayer2DDTI::Create: No source" );
  
  // Get some values from the MRI.
  float RASBounds[6];
  this->GetRASBounds( RASBounds );

  mWorldSize[0] = RASBounds[1] - RASBounds[0];
  mWorldSize[1] = RASBounds[3] - RASBounds[2];
  mWorldSize[2] = RASBounds[5] - RASBounds[4];

  mWorldCenter[0] = mDTIProperties->GetFAVolumeSource()->GetRASCenterX();
  mWorldCenter[1] = mDTIProperties->GetFAVolumeSource()->GetRASCenterY();
  mWorldCenter[2] = mDTIProperties->GetFAVolumeSource()->GetRASCenterZ();

  // 
  // Merged source
  //
  vtkImageAppendComponents* source = mDTIProperties->GetMergedSource();
  
  //
  // Transform
  //
  mVolumeToRAS = vtkImageReslice::New();
  mVolumeToRAS->InterpolateOn();

  // Nearest neighbor interpolation copies the vectors, so only linear
  // and cubic interpolations are allowed.
  mVolumeToRAS->
    SetInterpolationMode( mDTIProperties->GetTensorInterpolationType() );

  mVolumeToRAS->SetInputConnection( source->GetOutputPort() );
  mVolumeToRAS->SetOutputDimensionality( 3 );

  // This rotates the volume to the proper orientation. From
  // ImageReslice: "applying a transform to the resampling grid (which
  // lies in the output coordinate system) is equivalent to applying
  // the inverse of that transform to the input volume."
  double* rtv = mDTIProperties->GetFAVolumeSource()->GetRASToVoxelMatrix();
    
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

  mVolumeToRAS->SetResliceTransform( transform );
  mVolumeToRAS->BorderOff();
  transform->Delete();

  // This sets our output extent.
  mVolumeToRAS->SetOutputExtent( (int)RASBounds[0], (int)RASBounds[1],
                                (int)RASBounds[2], (int)RASBounds[3],
                                (int)RASBounds[4], (int)RASBounds[5] );
      
  // set up the transform for reorienting the tensors
  double *vtr = mDTIProperties->GetFAVolumeSource()->GetVoxelToRASMatrix();
  
  vtkMatrix4x4* voxelToRAS = vtkMatrix4x4::New();
  voxelToRAS->SetElement( 0, 0, vtr[0] * mDTIProperties->GetFAVolumeSource()->GetPixelSizeX() );
  voxelToRAS->SetElement( 0, 1, vtr[1] );
  voxelToRAS->SetElement( 0, 2, vtr[2] );
  voxelToRAS->SetElement( 0, 3, 0 );
  voxelToRAS->SetElement( 1, 0, vtr[4] );
  voxelToRAS->SetElement( 1, 1, vtr[5] * mDTIProperties->GetFAVolumeSource()->GetPixelSizeY() );
  voxelToRAS->SetElement( 1, 2, vtr[6] );
  voxelToRAS->SetElement( 1, 3, 0 );
  voxelToRAS->SetElement( 2, 0, vtr[8] );
  voxelToRAS->SetElement( 2, 1, vtr[9] );
  voxelToRAS->SetElement( 2, 2, vtr[10] * mDTIProperties->GetFAVolumeSource()->GetPixelSizeZ() );
  voxelToRAS->SetElement( 2, 3, 0 );
  voxelToRAS->SetElement( 3, 0, 0 );
  voxelToRAS->SetElement( 3, 1, 0 );
  voxelToRAS->SetElement( 3, 2, 0 );
  voxelToRAS->SetElement( 3, 3, 1 );
  
  vtkTransform *voxelTransform = vtkTransform::New();
  voxelTransform->PostMultiply();
  voxelTransform->SetMatrix( voxelToRAS );  
  voxelToRAS->Delete();

  //
  // The reslice object just takes a slice out of the volume.
  //
  if ( !mVolumeToRASSlice )
    mVolumeToRASSlice = vtkImageReslice::New();
  mVolumeToRASSlice->SetInputConnection( mVolumeToRAS->GetOutputPort() );
  mVolumeToRASSlice->BorderOff();
  
  // This sets us to extract slices.
  mVolumeToRASSlice->SetOutputDimensionality( 2 );
  
  // This will change depending what orienation we're in.
  mVolumeToRASSlice->SetResliceAxesDirectionCosines( 
    1, 0, 0,
    0, 1, 0,
    0, 0, 1 );

  // This will change to select a different slice.
  mVolumeToRASSlice->SetResliceAxesOrigin( 0, 0, 0 );

  //
  // Flip over the x axis (left/right). This get us into neurological
  // view.
  //
  vtkImageFlip* imageFlip = vtkImageFlip::New();
  imageFlip->SetInputConnection( mVolumeToRASSlice->GetOutputPort() );
  imageFlip->SetFilteredAxis( 0 ); // x axis
  
  //
  // Glyph
  //
  mGlyph = vtkFDTensorGlyph::New();
  mGlyph->SetInputConnection( imageFlip->GetOutputPort() );
  imageFlip->Delete();

  // this reorients each tensor because when we resample the volume, the vector
  // data isn't reoriented
  mGlyph->SetVoxelToMeasurementFrameTransform( voxelTransform );  
  voxelTransform->Delete();
    
  //
  // Plane mapper transform.
  //
  if ( !mPlaneTransform )
    mPlaneTransform = vtkTransform::New();

  vtkTransformPolyDataFilter* planePDF = vtkTransformPolyDataFilter::New();
  planePDF->SetInput( mGlyph->GetOutput() );
  planePDF->SetTransform( mPlaneTransform );
      
  //
  // Mapper
  //
  mMapper = vtkPolyDataMapper::New();
  mMapper->UseLookupTableScalarRangeOn();
  mMapper->SetInput( planePDF->GetOutput() );
  mMapper->SetLookupTable( (vtkScalarsToColors*)mGlyph->GetColorTable() );
  
  planePDF->Delete();
  
  //
  // Main Actor
  //
  mActor = vtkLODActor::New();
  mActor->SetMapper( mMapper );
  mActor->GetProperty()->SetAmbient( .3 );
  mActor->GetProperty()->SetDiffuse( .2 );
  mActor->GetProperty()->SetSpecular( .2 );

  // Set the size and number of points for the LOD when it goes into
  // point mode.
  mActor->GetProperty()->SetPointSize( 5.0 );
  mActor->SetNumberOfCloudPoints( 800 );
  
  // 
  // Edge Mapper
  //
  vtkPolyDataMapper* edgeMapper = vtkPolyDataMapper::New();
  edgeMapper->UseLookupTableScalarRangeOn();
  edgeMapper->SetLookupTable( (vtkScalarsToColors*)mGlyph->GetColorTable() );
  edgeMapper->SetInput( planePDF->GetOutput() );

  //
  // Edge Actor
  //
  mEdgeActor = vtkLODActor::New();
  mEdgeActor->SetMapper( edgeMapper );
  mEdgeActor->GetProperty()->SetRepresentationToWireframe();
  mEdgeActor->GetProperty()->SetColor( 0, 0, 0 );
  mEdgeActor->GetProperty()->SetLineWidth( 1.4 );

  // Initial visibility.
  if( !mDTIProperties->GetRenderEdges() )
    mEdgeActor->SetVisibility( false );

  this->AddProp( mActor );
  this->AddProp( mEdgeActor );
  
  // Set up for the initial viewing state.
  this->UpdateOpacity();
  this->UpdateDetail();
  this->UpdateGlyphScaling();
  this->Update2DInfo();
}

void
vtkKWScubaLayer2DDTI::AddControls ( vtkKWWidget* iPanel ) {
  
}

void
vtkKWScubaLayer2DDTI::RemoveControls () {
  
}

void
vtkKWScubaLayer2DDTI::DoListenToMessage ( string isMessage, void* iData ) {

  if( isMessage == "OpacityChanged" ) {
    this->UpdateOpacity();
  } else if( isMessage == "TensorDetailChanged" ) {
    this->UpdateDetail();
  } else if( isMessage == "TensorScalingChanged" ) {
    this->UpdateGlyphScaling();
  } else if( isMessage == "RenderEdgesChanged" ) {
    this->UpdateEdges();
  } else if( isMessage == "Layer2DInfoChanged" ) {
    this->Update2DInfo();
  } else if( isMessage == "TensorInterpolationChanged" ) {

    if( NULL != mDTIProperties && NULL != mVolumeToRAS ) {
      
      mVolumeToRAS->SetInterpolationMode( 
        mDTIProperties->GetTensorInterpolationType() );

      this->PipelineChanged();
    }

  }
}

void
vtkKWScubaLayer2DDTI::GetRASBounds ( float ioBounds[6] ) const {

  if( mDTIProperties && mDTIProperties->GetFAVolumeSource() )
    mDTIProperties->GetFAVolumeSource()->GetRASBounds( ioBounds );
  else
    vtkKWScubaLayer::GetRASBounds( ioBounds );
}

void
vtkKWScubaLayer2DDTI::Get2DRASZIncrementHint ( float ioHint[3]) const {

  ioHint[0] = 1.0;
  ioHint[1] = 1.0;
  ioHint[2] = 1.0;
}

void
vtkKWScubaLayer2DDTI::GetInfoItems ( float iRAS[3],
                                      list<ScubaInfoItem>& ilInfo ) const {
                                        
  if( mDTIProperties ) {

    // get the index of the RAS our FA volume
    vtkFSVolumeSource* faSource = mDTIProperties->GetFAVolumeSource();
    int index[ 3 ];
    faSource->ConvertRASToIndex( iRAS[0], iRAS[1], iRAS[2],
      index[0], index[1], index[2] );
      
    // get the fa at this index
    const float fa = faSource->GetValueAtIndex( index[0], index[1], index[2] );
    
    char infoValue[1024];
    snprintf( infoValue, sizeof(infoValue), "%f", fa );

    ScubaInfoItem faInfo;
    faInfo.Clear();
    faInfo.SetLabel( "FA" );
    faInfo.SetValue( infoValue );
    faInfo.SetShortenHint( false );

    // Return it.
    ilInfo.push_back( faInfo );
 
    // return the eigenvalues
    vtkFSVolumeSource* eigenValuesSource = mDTIProperties->GetEigenValueVolumeSource();
      
    // each location will have three eigenvalues
    float eigenValues[ 3 ];
    for( int cValue=0; cValue<3; cValue++ ) {
      eigenValues[ cValue ] = eigenValuesSource->GetValueAtIndex( index[0], index[1], index[2], cValue );
    }
    
    snprintf( infoValue, sizeof(infoValue), "%f %f %f", eigenValues[ 0 ], eigenValues[ 1 ], eigenValues[ 2 ] );
    ScubaInfoItem eigenValueInfo;
    eigenValueInfo.Clear();
    eigenValueInfo.SetLabel( "Eigen Values" );
    eigenValueInfo.SetValue( infoValue );
    eigenValueInfo.SetShortenHint( false );

    // Return the eigenvalues
    ilInfo.push_back( eigenValueInfo );  
    
    for( int n=0; n<3; n++ ) {
      
      const int currentEigenVectorId = n + 1;
      
      // get the eigenvectors
      vtkFSVolumeSource* eigenVectorSource = mDTIProperties->GetEigenVectorVolumeSource( currentEigenVectorId );
      float eigenVector[ 3 ];
  
      for( int cValue=0; cValue<3; cValue++ ) {
        eigenVector[ cValue ] = eigenVectorSource->GetValueAtIndex( index[0], index[1], index[2], cValue );
      }
      
      snprintf( infoValue, sizeof(infoValue), "%f %f %f", eigenVector[ 0 ], eigenVector[ 1 ], eigenVector[ 2 ] );
      ScubaInfoItem eigenVectorInfo;
      eigenVectorInfo.Clear();
      
      char label[1024];
      snprintf( label, sizeof(label), "Eigen Vector %i", currentEigenVectorId );
      eigenVectorInfo.SetLabel( label );
      eigenVectorInfo.SetValue( infoValue );
      eigenVectorInfo.SetShortenHint( false );
  
      // Return the eigenvector
      ilInfo.push_back( eigenVectorInfo );  
    }

  
    
  } else {

    ScubaInfoItem info;
    info.Clear();
    info.SetLabel( "" );
    info.SetValue( "" );
    info.SetShortenHint( false );
  
    // Return it.
    ilInfo.push_back( info );
    
  }
}

void
vtkKWScubaLayer2DDTI::UpdateOpacity () {
  if ( mActor )
    if ( mActor->GetProperty() ) {
      mActor->GetProperty()->SetOpacity( mProperties->GetOpacity() );
      this->PipelineChanged();
    }
}

void
vtkKWScubaLayer2DDTI::UpdateDetail () {

  if( mVolumeToRAS ) {
    
    const int shrinkage = this->GetCurrentShrinkageValue();
    
    // this changes the pixel/tensor sampling
    mVolumeToRAS->SetOutputSpacing(
      mDTIProperties->GetFAVolumeSource()->GetPixelSizeX() * shrinkage,
      mDTIProperties->GetFAVolumeSource()->GetPixelSizeY() * shrinkage,
      mDTIProperties->GetFAVolumeSource()->GetPixelSizeZ() * shrinkage );
        
    this->Update2DInfo();
  }
}

void
vtkKWScubaLayer2DDTI::UpdateGlyphScaling () {

  bool bUniform = false;
  bool bFA = false;
  
  switch( mDTIProperties->GetTensorScaling() ) {
  case ScubaCollectionPropertiesDTI::Uniform:
    bUniform = true;
    break;
  case ScubaCollectionPropertiesDTI::FA:
    bFA = true;
    break;
  }
    
  mGlyph->SetUniformScaling( bUniform );
  mGlyph->SetFAScaling( bFA );
  
  this->PipelineChanged();
}

void
vtkKWScubaLayer2DDTI::Update2DInfo () {
  
  if( NULL == mProperties )
    return;
    
  // go through the plane
  const int inPlane = mViewProperties->Get2DInPlane();

  // current slice
  const float rasZ = mViewProperties->Get2DRASZ();

  // our slices will be rotated, so we need to reorient the tensors again
  vtkTransform *sliceTransform = vtkTransform::New();
  sliceTransform->Identity();
  
  // new orientations for each plane
  if( inPlane == 0 ) {
    sliceTransform->RotateY( 90 );
    sliceTransform->RotateX( 90 );
    sliceTransform->RotateZ( 180 );
  } else if( inPlane == 1 ) {
    sliceTransform->RotateY( 180 );
    sliceTransform->RotateX( 90 );
  } else if( inPlane == 2 ) {
    sliceTransform->RotateY( 180 );
    sliceTransform->RotateX( 180 );
  }
  mGlyph->SetVoxelToSliceTransform( sliceTransform );
  sliceTransform->Delete();

  //
  // Now we'll set up our three planes with the right transforms.
  //
  double translation[ 3 ] = { mWorldCenter[ 0 ], mWorldCenter[ 1 ], mWorldCenter[ 2 ] };
  translation[ inPlane ] = 0;
  
  mPlaneTransform->Identity();
  mPlaneTransform->Translate( translation );
 
  // the direction of our reslice.  We'll want to change them for each
  // plane.  The cosine direction is the vector perpendicular to the plane
  double directionCosines[9] = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
  
  // rotate our slices to get them in the right orientation
  if( inPlane == 0 ) {
    
    mPlaneTransform->RotateX( 90 );
    mPlaneTransform->RotateY( 90 );
    
    // Putting negatives in the reslice axes cosines will flip the
    // image on that axis.
    directionCosines[ 1 ] = -1;
    directionCosines[ 5 ] = 1;
    directionCosines[ 6 ] = 1;
    
    
  } else if( inPlane == 1 ) {

    mPlaneTransform->RotateX( 90 );
    mPlaneTransform->RotateY( 180 );
    
    directionCosines[ 0 ] = 1;
    directionCosines[ 5 ] = 1;
    directionCosines[ 7 ] = 1;
    
  } else if( inPlane == 2 ) {

    mPlaneTransform->RotateY( 180 );

    directionCosines[ 0 ] = 1;
    directionCosines[ 4 ] = 1;
    directionCosines[ 8 ] = 1;
    
  }

  mVolumeToRASSlice->SetResliceAxesDirectionCosines( directionCosines );
  
  // set the origin for each plane
  double axesOrigin[3] = { 0, 0, 0 };      
  axesOrigin[ inPlane ] = static_cast< int >( rasZ - mWorldCenter[ inPlane ] );
  mVolumeToRASSlice->SetResliceAxesOrigin( axesOrigin );

  this->PipelineChanged();
  
}

int
vtkKWScubaLayer2DDTI::GetCurrentShrinkageValue () {
  
  int shrinkage = 0;

  switch( mDTIProperties->GetTensorDetail() ) {
  case ScubaCollectionPropertiesDTI::Least:
    shrinkage = 4;
    break;
  case ScubaCollectionPropertiesDTI::Less:
    shrinkage = 2;
    break;
  case ScubaCollectionPropertiesDTI::Normal:
    shrinkage = 1;
    break;
  }
  
  return shrinkage;
}

void
vtkKWScubaLayer2DDTI::UpdateEdges() {
  
  const bool isEdgeVisible = mDTIProperties->GetRenderEdges();
  mEdgeActor->SetVisibility( isEdgeVisible );    

  this->PipelineChanged();
}

