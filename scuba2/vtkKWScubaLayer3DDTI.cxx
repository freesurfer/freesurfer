/**
 * @file  vtkKWScubaLayer3DDTI.cxx
 * @brief Displays 3D glyphs
 *
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:40 $
 *    $Revision: 1.10 $
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

#include "vtkKWScubaLayer3DDTI.h"

#include <string>
#include <stdexcept>

#include "ScubaCollectionProperties.h"
#include "ScubaCollectionPropertiesDTI.h"
#include "ScubaInfoItem.h"
#include "vtkLODActor.h"
#include "vtkFDTensorGlyph.h"
#include "vtkImageAppendComponents.h"
#include "vtkMatrix4x4.h"
#include "vtkObjectFactory.h"
#include "vtkPolyDataMapper.h"
#include "vtkProperty.h"
#include "vtkTransform.h"
#include "vtkFSVolumeSource.h"
#include "vtkImageReslice.h"
#include "vtkImageShrink3D.h"
#include "vtkTransformPolyDataFilter.h"
#include "vtkImageFlip.h"

#include "vtkKWCheckButton.h"
#include "vtkKWRadioButtonSet.h"
#include "vtkKWRadioButton.h"


using namespace std;

vtkStandardNewMacro( vtkKWScubaLayer3DDTI );
vtkCxxRevisionMacro( vtkKWScubaLayer3DDTI, "$Revision: 1.10 $" );

vtkKWScubaLayer3DDTI::vtkKWScubaLayer3DDTI () :
  mDTIProperties( NULL ),
  mVolumeToRAS( NULL ),
  mIsPlaneXVisbleButton( NULL ),
  mIsPlaneYVisbleButton( NULL ),
  mIsPlaneZVisbleButton( NULL ),
  mIsShowingAllButton( NULL ),
  mIsPlaneXVisible( true ),
  mIsPlaneYVisible( true ),
  mIsPlaneZVisible( true ),
  mIsShowingAll( false ) {
    
  for( int n = 0; n < 3; n++ ) {
    mGlyphs[ n ] = NULL;
    mVolumeToRASSlice[ n ] = NULL;
    mSliceTransform[ n ] = NULL;
    mPlaneMappers[ n ] = NULL;
    mPlaneActors[ n ] = NULL;
    mGlyphEdgeActors[ n ] = NULL;
  }
    
}

vtkKWScubaLayer3DDTI::~vtkKWScubaLayer3DDTI () {
  
  if ( mVolumeToRAS )
    mVolumeToRAS->Delete();
    
  for( int n=0; n<3; n++ ) {
  
    if ( mGlyphs[ n ] )
      mGlyphs[ n ]->Delete();
      
    if( mVolumeToRASSlice[ n ] )
      mVolumeToRASSlice[ n ]->Delete();

    if( mSliceTransform[ n ] )
      mSliceTransform[ n ]->Delete();

    if ( mPlaneMappers[ n ] )
      mPlaneMappers[ n ]->Delete();

    if ( mPlaneActors[ n ] )
      mPlaneActors[ n ]->Delete();

    if ( mGlyphEdgeActors[ n ] )
      mGlyphEdgeActors[ n ]->Delete();
      
  }
  
}

void
vtkKWScubaLayer3DDTI::SetDTIProperties ( ScubaCollectionPropertiesDTI const* iProperties ) {
  mDTIProperties = iProperties;
}


void
vtkKWScubaLayer3DDTI::Create () {

  // Bail if we don't have our source yet.
  if( NULL == mDTIProperties )
    throw runtime_error( "vtkKWScubaLayer3DDTI::Create: No source" );

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

  // set up the transform for reorienting the tensors
  double vtr[ 16 ];
  mDTIProperties->GetFAVolumeSource()->GetUnscaledVoxelToRASMatrix( vtr );

  vtkMatrix4x4* voxelToRAS = vtkMatrix4x4::New();
  voxelToRAS->SetElement( 0, 0, vtr[0] );
  voxelToRAS->SetElement( 0, 1, vtr[1] );
  voxelToRAS->SetElement( 0, 2, vtr[2] );
  voxelToRAS->SetElement( 0, 3, 0 );

  voxelToRAS->SetElement( 1, 0, vtr[4] );
  voxelToRAS->SetElement( 1, 1, vtr[5] );
  voxelToRAS->SetElement( 1, 2, vtr[6] );
  voxelToRAS->SetElement( 1, 3, 0 );

  voxelToRAS->SetElement( 2, 0, vtr[8] );
  voxelToRAS->SetElement( 2, 1, vtr[9] );
  voxelToRAS->SetElement( 2, 2, vtr[10] );
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
  // This transforms the voxel space source into RAS space.
  //
  mVolumeToRAS = vtkImageReslice::New();
  mVolumeToRAS->InterpolateOn();
  
  // the default interpolation (nearest neighbor) copies the vectors, so only
  // linear or cubic interoplations are allowed
  mVolumeToRAS->SetInterpolationMode( 
    mDTIProperties->GetTensorInterpolationType() );

  mVolumeToRAS->SetInputConnection( source->GetOutputPort() );
  mVolumeToRAS->SetOutputDimensionality( 3 );
   
  double *rtv = mDTIProperties->GetFAVolumeSource()->GetRASToVoxelMatrix();
  
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
  
  // three times for three cut planes
  for( int n = 0; n<3; n++ ) {  
      
    //
    // The reslice object just takes a slice out of the volume.
    //
    if ( !mVolumeToRASSlice[n] )
      mVolumeToRASSlice[n] = vtkImageReslice::New();
    mVolumeToRASSlice[n]->SetInputConnection( mVolumeToRAS->GetOutputPort());
    mVolumeToRASSlice[n]->BorderOff();
    
    // This sets us to extract slices.
    mVolumeToRASSlice[n]->SetOutputDimensionality( 2 );
    
    // This will change depending what orienation we're in.
    mVolumeToRASSlice[n]->SetResliceAxesDirectionCosines( 1, 0, 0,
             0, 1, 0,
             0, 0, 1 );
    
    // This will change to select a different slice.
    mVolumeToRASSlice[n]->SetResliceAxesOrigin( 0, 0, 0 );
  
    //
    // Flip over the x axis (left/right). This get us into neurological
    // view.
    //
    vtkImageFlip* imageFlip = vtkImageFlip::New();
    imageFlip->SetInputConnection( mVolumeToRASSlice[n]->GetOutputPort() );
    imageFlip->SetFilteredAxis( 0 ); // x axis

    mGlyphs[ n ] = vtkFDTensorGlyph::New();
    mGlyphs[ n ]->SetInputConnection( imageFlip->GetOutputPort() );
    
    // this reorients each tensor because when we resample the volume, the vector
    // data isn't reoriented
    mGlyphs[ n ]->SetVoxelToMeasurementFrameTransform( voxelTransform );

    imageFlip->Delete();
    
    // our slices will be rotated, so we need to reorient the tensors again
    vtkTransform *sliceTransform = vtkTransform::New();
    sliceTransform->Identity();
    
    // new orientations for each plane
    if( n == 0 ) {
      sliceTransform->RotateY( 90 );
      sliceTransform->RotateX( 90 );
      sliceTransform->RotateZ( 180 );
    } else if( n == 1 ) {
      sliceTransform->RotateY( 180 );
      sliceTransform->RotateX( 90 );
    } else if( n == 2 ) {
      sliceTransform->RotateY( 180 );
      sliceTransform->RotateX( 180 );
    }
    mGlyphs[ n ]->SetVoxelToSliceTransform( sliceTransform );
    sliceTransform->Delete();    

    // Plane mapper transform for moving the slice plane
    if ( !mSliceTransform[n] )
      mSliceTransform[n] = vtkTransform::New();
    
    //
    // Poly data from plane and plane transform.
    //
    vtkTransformPolyDataFilter* planePDF = vtkTransformPolyDataFilter::New();
    planePDF->SetInput( mGlyphs[ n ]->GetOutput() );
    planePDF->SetTransform( mSliceTransform[n] );

    mPlaneMappers[ n ] = vtkPolyDataMapper::New();  
    mPlaneMappers[ n ]->UseLookupTableScalarRangeOn();
    mPlaneMappers[ n ]->SetLookupTable( (vtkScalarsToColors*)mGlyphs[ n ]->GetColorTable() );
    mPlaneMappers[ n ]->SetInputConnection( planePDF->GetOutputPort() );
    
    mPlaneActors[ n ] = vtkLODActor::New();
    mPlaneActors[ n ]->SetMapper( mPlaneMappers[ n ] );
    mPlaneActors[ n ]->GetProperty()->SetAmbient( .3 );
    mPlaneActors[ n ]->GetProperty()->SetDiffuse( .2 );
    mPlaneActors[ n ]->GetProperty()->SetSpecular( .2 );
    
    // this set the size and number of points for the level of detail when it
    // goes in point mode
    mPlaneActors[ n ]->GetProperty()->SetPointSize( 5.0 );
    mPlaneActors[ n ]->SetNumberOfCloudPoints( 800 );
    
    // create the edges
    vtkPolyDataMapper* edgePlaneMappers = vtkPolyDataMapper::New();  
    edgePlaneMappers->UseLookupTableScalarRangeOn();
    edgePlaneMappers->SetLookupTable( (vtkScalarsToColors*)mGlyphs[ n ]->GetColorTable() );    
    edgePlaneMappers->SetInputConnection( planePDF->GetOutputPort() );

    planePDF->Delete();

    mGlyphEdgeActors[ n ] = vtkLODActor::New();
    mGlyphEdgeActors[ n ]->SetMapper( edgePlaneMappers );

    vtkProperty* property = mGlyphEdgeActors[ n ]->GetProperty();
    property->SetRepresentationToWireframe();
    property->SetColor( 0, 0, 0 );
    property->SetLineWidth( 1.4 );
    
    if( !mDTIProperties->GetRenderEdges() ) {
      mGlyphEdgeActors[ n ]->SetVisibility( false );
    }

    this->AddProp( mPlaneActors[ n ] );
    this->AddProp( mGlyphEdgeActors[ n ] );
  }

  voxelTransform->Delete();
  
  this->UpdateDetail();
  this->UpdateGlyphScaling();
  this->UpdatePlanes();
    
}

void
vtkKWScubaLayer3DDTI::AddControls ( vtkKWWidget* iPanel ) {

  // create a check boxes for toggling the display of the tensor slices
  mIsPlaneXVisbleButton = vtkKWCheckButton::New();
  mIsPlaneXVisbleButton->SetParent( iPanel );
  mIsPlaneXVisbleButton->Create();
  mIsPlaneXVisbleButton->SetAnchorToWest();
  mIsPlaneXVisbleButton->SetText( "Show X RAS Tensors" );
  mIsPlaneXVisbleButton->SetSelectedState( mIsPlaneXVisible );  
  // set the callback
  mIsPlaneXVisbleButton->SetCommand( this, "SetIsPlaneXVisible" );
  mIsPlaneXVisbleButton->SetEnabled( !mIsShowingAll );  
  
  mIsPlaneYVisbleButton = vtkKWCheckButton::New();
  mIsPlaneYVisbleButton->SetParent( iPanel );
  mIsPlaneYVisbleButton->Create();
  mIsPlaneYVisbleButton->SetAnchorToWest();
  mIsPlaneYVisbleButton->SetText( "Show Y RAS Tensors" );
  mIsPlaneYVisbleButton->SetSelectedState( mIsPlaneYVisible );  
  // set the callback
  mIsPlaneYVisbleButton->SetCommand( this, "SetIsPlaneYVisible" );
  mIsPlaneYVisbleButton->SetEnabled( !mIsShowingAll );  

  mIsPlaneZVisbleButton = vtkKWCheckButton::New();
  mIsPlaneZVisbleButton->SetParent( iPanel );
  mIsPlaneZVisbleButton->Create();
  mIsPlaneZVisbleButton->SetAnchorToWest();
  mIsPlaneZVisbleButton->SetText( "Show Z RAS Tensors" );
  mIsPlaneZVisbleButton->SetSelectedState( mIsPlaneZVisible );  
  // set the callback
  mIsPlaneZVisbleButton->SetCommand( this, "SetIsPlaneZVisible" );
  mIsPlaneZVisbleButton->SetEnabled( !mIsShowingAll );  
  
  mIsShowingAllButton = vtkKWCheckButton::New();
  mIsShowingAllButton->SetParent( iPanel );
  mIsShowingAllButton->Create();
  mIsShowingAllButton->SetAnchorToWest();
  mIsShowingAllButton->SetText( "Show all tensors" );
  mIsShowingAllButton->SetSelectedState( mIsShowingAll );  
  // set the callback
  mIsShowingAllButton->SetCommand( this, "SetIsShowingAll" );
  
  // TODO: disabling the show all button because nothing is implemented for it now
  mIsShowingAllButton->SetEnabled( false );  
    
  this->Script( "pack %s %s %s %s -side top -fill x -anchor nw",
                mIsPlaneXVisbleButton->GetWidgetName(), 
                mIsPlaneYVisbleButton->GetWidgetName(),
                mIsPlaneZVisbleButton->GetWidgetName(),
                mIsShowingAllButton->GetWidgetName() );

}

void
vtkKWScubaLayer3DDTI::RemoveControls () {
  mIsPlaneXVisbleButton->Delete();
  mIsPlaneXVisbleButton = NULL;
}

void
vtkKWScubaLayer3DDTI::DoListenToMessage ( string isMessage, void* iData ) {

  if( isMessage == "OpacityChanged" ) {
    
    for( int n=0; n<3; n++ ) {

      if ( mPlaneActors[ n ] )
        if ( mPlaneActors[ n ]->GetProperty() ) {
          mPlaneActors[ n ]->GetProperty()->SetOpacity( mProperties->GetOpacity() );
          this->PipelineChanged();
        }
    }

  } else if( isMessage == "RenderEdgesChanged" ) {

    this->UpdateEdges();

  } else if( isMessage == "Layer3DInfoChanged" ) {

    this->UpdatePlanes();

  } else if( isMessage == "FastModeChanged" ) {

    if ( mViewProperties->GetFastMode() ) {

      // we probably don't need to do anything with less detail
      if( mDTIProperties->GetTensorDetail() != ScubaCollectionPropertiesDTI::Least ) {

        for( int n = 0; n<3; n++ ) {          
          mPlaneActors[ n ]->GetProperty()->SetRepresentationToPoints ();
        }
        
      }
      
    } else {

      for( int n = 0; n<3; n++ ) {  
        if( VTK_POINTS == mPlaneActors[ n ]->GetProperty()->GetRepresentation() ) {
          mPlaneActors[ n ]->GetProperty()->SetRepresentationToSurface ();
        }
      }
      
    }
    
  } else if( isMessage == "TensorInterpolationChanged" ) {

    if( NULL != mDTIProperties && NULL != mVolumeToRAS ) {
      mVolumeToRAS->SetInterpolationMode( 
        mDTIProperties->GetTensorInterpolationType() );
    }
      
    this->PipelineChanged();
    
  } else if( isMessage == "TensorScalingChanged" ) {

    this->UpdateGlyphScaling();

  } else if( isMessage == "TensorDetailChanged" ) {

    this->UpdateDetail();

  }

  
}

void
vtkKWScubaLayer3DDTI::GetRASBounds ( float ioBounds[6] ) const {
  if ( mDTIProperties ) {
    mDTIProperties->GetFAVolumeSource()->GetRASBounds( ioBounds );
  } else {
    for ( int nBound = 0; nBound < 6; nBound++ )
      ioBounds[nBound] = 0;
  }
}

void
vtkKWScubaLayer3DDTI::Get3DRASZIncrementHint ( float ioHint[3]) const {

  ioHint[0] = 1.0;
  ioHint[1] = 1.0;
  ioHint[2] = 1.0;
}

void
vtkKWScubaLayer3DDTI::GetInfoItems ( float iRAS[3],
                                      list<ScubaInfoItem>& ilInfo ) const {

  ScubaInfoItem info;

  if ( mDTIProperties ) {

    // get the index of the RAS our FA volume
    vtkFSVolumeSource* faSource = mDTIProperties->GetFAVolumeSource();
    int index[ 3 ];
    faSource->ConvertRASToIndex( iRAS[0], iRAS[1], iRAS[2],
      index[0], index[1], index[2] );
      
    // get the fa at this index
    const float fa = faSource->GetValueAtIndex( index[0], index[1], index[2] );
    
    char infoValue[1024];
    snprintf( infoValue, sizeof(infoValue), "%f", fa );

    info.Clear();
    info.SetLabel( "FA" );
    info.SetValue( infoValue );
    info.SetShortenHint( false );
  
    // Return the fa
    ilInfo.push_back( info );  
    
    // get the eigenvalues
    //    vtkFSVolumeSource* eigenValuesSource = mDTIProperties->GetEigenValueVolumeSource();
//    eigenValuesSource->ConvertRASToIndex( iRAS[0], iRAS[1], iRAS[2],
//      index[0], index[1], index[2] );
//      
//    // each location will have three eigenvalues
//    float eigenValues[ 3 ];
//    for( int cValue=0; cValue<3; cValue++ ) {
//      cerr << "  getting value: " << cValue << std::endl;
//      eigenValues[ cValue ] = eigenValuesSource->GetValueAtIndex( index[0], index[1], index[2], cValue );
//    }
//    snprintf( infoValue, sizeof(infoValue), "%f %f %f", eigenValues[ 0 ], eigenValues[ 1 ], eigenValues[ 2 ] );
//    info.Clear();
//    info.SetLabel( "Eigen Values" );
//    info.SetValue( infoValue );
//    info.SetShortenHint( false );
//
//    // Return the eigenvalues
//    ilInfo.push_back( info );  
    
      
//    vtkFSVolumeSource* eigenVector1Source = mDTIProperties->GetEigenVector1VolumeSource();
//    vtkFSVolumeSource* eigenVector2Source = mDTIProperties->GetEigenVector2VolumeSource();
//    vtkFSVolumeSource* eigenVector3Source = mDTIProperties->GetEigenVector3VolumeSource();
      
  } else {
    info.Clear();
    info.SetLabel( "" );
    info.SetValue( "" );
    info.SetShortenHint( false );
  
    // Return it.
    ilInfo.push_back( info );  
  }  
  
}

void 
vtkKWScubaLayer3DDTI::SetIsPlaneXVisible ( int ibIsVisible ) {
  if( ibIsVisible != mIsPlaneXVisible ) {
    mIsPlaneXVisible = ibIsVisible;
    mPlaneActors[ X_PLANE ]->SetVisibility( mIsPlaneXVisible );
    this->PipelineChanged();
  }
  this->UpdateEdges();
}

void 
vtkKWScubaLayer3DDTI::SetIsPlaneYVisible ( int ibIsVisible ) {
  if( ibIsVisible != mIsPlaneYVisible ) {
    mIsPlaneYVisible = ibIsVisible;
    mPlaneActors[ Y_PLANE ]->SetVisibility( mIsPlaneYVisible );
    this->PipelineChanged();
  }
  this->UpdateEdges();
}

void 
vtkKWScubaLayer3DDTI::SetIsPlaneZVisible ( int ibIsVisible ) {
  if( ibIsVisible != mIsPlaneZVisible ) {
    mIsPlaneZVisible = ibIsVisible;
    mPlaneActors[ Z_PLANE ]->SetVisibility( mIsPlaneZVisible );
    this->PipelineChanged();
  }
  this->UpdateEdges();
}

void 
vtkKWScubaLayer3DDTI::SetIsShowingAll( int ibIsShowingAll ) {

  if( mIsShowingAll != ibIsShowingAll ) {
    
    mIsShowingAll = ibIsShowingAll;
    
    mIsPlaneXVisbleButton->SetEnabled( !mIsShowingAll );
    mIsPlaneYVisbleButton->SetEnabled( !mIsShowingAll );
    mIsPlaneZVisbleButton->SetEnabled( !mIsShowingAll );
    
    UpdatePlanes();
    
  }
}

void 
vtkKWScubaLayer3DDTI::UpdatePlanes () {
  this->UpdatePlanes( false );
}

void 
vtkKWScubaLayer3DDTI::UpdatePlanes ( bool hasDetailChanged ) {
    
  if( mIsShowingAll ) {
    
    // TODO: implement something here?

  } else {
    
    const int ras[ 3 ] = {
      static_cast< int >( mViewProperties->Get3DRASX() ),
      static_cast< int >( mViewProperties->Get3DRASY() ),
      static_cast< int >( mViewProperties->Get3DRASZ() )
    };
    
    // go through each plane
    for( int n=0; n<3; n++ ) {

      //
      // Now we'll set up our three planes with the right transforms.
      //
      double translation[ 3 ] = { 0, 0, 0 };
      translation[ n ] = ras[ n ] - mWorldCenter[ n ];
      
      mSliceTransform[ n ]->Identity();
      mSliceTransform[ n ]->Translate( translation );
     
      // the direction of our reslice.  We'll want to change them for each
      // plane.  The cosine direction is the vector perpendicular to the plane
      double directionCosines[9] = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
      
      // rotate our slices to get them in the right orientation
      if( n == 0 ) {
        
        mSliceTransform[n]->RotateX( 90 );
        mSliceTransform[n]->RotateY( 90 );
        
        // Putting negatives in the reslice axes cosines will flip the
        // image on that axis.
        directionCosines[ 1 ] = -1;
        directionCosines[ 5 ] = 1;
        directionCosines[ 6 ] = 1;
        
        
      } else if( n == 1 ) {

        mSliceTransform[n]->RotateX( 90 );
        mSliceTransform[n]->RotateY( 180 );
        
        directionCosines[ 0 ] = 1;
        directionCosines[ 5 ] = 1;
        directionCosines[ 7 ] = 1;
        
      } else if( n == 2 ) {

        mSliceTransform[n]->RotateY( 180 );

        directionCosines[ 0 ] = 1;
        directionCosines[ 4 ] = 1;
        directionCosines[ 8 ] = 1;
        
      }

      mVolumeToRASSlice[n]->SetResliceAxesDirectionCosines( directionCosines );
      
      // set the origin for each plane
      double axesOrigin[3] = { 0, 0, 0 };      
      
      // we should take the floor here, not round, truncate, or take the ceiling 
      axesOrigin[ n ] = static_cast< int >( floor( ras[n] - mWorldCenter[n] ) );
      mVolumeToRASSlice[n]->SetResliceAxesOrigin( axesOrigin );
            
    }
                      
  }
  
  this->PipelineChanged();
  
}

void
vtkKWScubaLayer3DDTI::UpdateGlyphScaling () {
  bool isUniform = false;
  bool isFA = false;
  
  switch( mDTIProperties->GetTensorScaling() ) {
  case ScubaCollectionPropertiesDTI::Uniform:
    isUniform = true;
    break;
  case ScubaCollectionPropertiesDTI::FA:
    isFA = true;
    break;
  }
    
  for( int n=0; n<3; n++ ) {
    mGlyphs[ n ]->SetUniformScaling( isUniform );
    mGlyphs[ n ]->SetFAScaling( isFA );
  }

  this->PipelineChanged();
}

void
vtkKWScubaLayer3DDTI::UpdateDetail () {
  
  if( mVolumeToRAS ) {

    const int shrinkage = this->GetCurrentShrinkage();
    
    // this changes the pixel/tensor sampling
    mVolumeToRAS->SetOutputSpacing(
      mDTIProperties->GetFAVolumeSource()->GetPixelSizeX() * shrinkage,
      mDTIProperties->GetFAVolumeSource()->GetPixelSizeY() * shrinkage,
      mDTIProperties->GetFAVolumeSource()->GetPixelSizeZ() * shrinkage );

    // if we update the detail, we also need to update the planes to be at the
    // positioned correctly
    this->UpdatePlanes( true );
    
  }
    
  this->PipelineChanged();
}

int
vtkKWScubaLayer3DDTI::GetCurrentShrinkage () {
  
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
vtkKWScubaLayer3DDTI::UpdateEdges() {

  const bool arePlanesVisible[] = 
    { mIsPlaneXVisible, mIsPlaneYVisible, mIsPlaneZVisible };
  
  for( int n=0; n<3; n++ ) {
    
    if( arePlanesVisible[ n ] ) {
      const bool isEdgeVisible = mDTIProperties->GetRenderEdges();
      mGlyphEdgeActors[ n ]->SetVisibility( isEdgeVisible );
    }
    
  }

  this->PipelineChanged();
}

