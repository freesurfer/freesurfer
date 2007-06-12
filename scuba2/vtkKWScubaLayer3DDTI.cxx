/**
 * @file  vtkKWScubaLayer3DDTI.cxx
 * @brief Displays 3D glyphs
 *
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: dsjen $
 *    $Date: 2007/06/12 15:46:43 $
 *    $Revision: 1.4 $
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

#include "vtkKWScubaLayer3DDTI.h"

#include <string>
#include <stdexcept>

#include "ScubaCollectionProperties.h"
#include "ScubaCollectionPropertiesDTI.h"
#include "ScubaInfoItem.h"
#include "vtkLODActor.h"
#include "vtkFDTensorGlyph.h"
#include "vtkImageAppendComponents.h"
#include "vtkImageClip.h"
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
vtkCxxRevisionMacro( vtkKWScubaLayer3DDTI, "$Revision: 1.4 $" );

vtkKWScubaLayer3DDTI::vtkKWScubaLayer3DDTI () :
  mDTIProperties( NULL ),
  mVolumeToRAS( NULL ),
  mReducedVolume( NULL ),
  mIsPlaneXVisbleButton( NULL ),
  mIsPlaneYVisbleButton( NULL ),
  mIsPlaneZVisbleButton( NULL ),
  mIsShowingAllButton( NULL ),
  mIsPlaneXVisible( true ),
  mIsPlaneYVisible( true ),
  mIsPlaneZVisible( true ),
  mIsShowingAll( false ) {
    
  for( int n = 0; n < 3; n++ ) {
    mClips[ n ] = NULL;
    mGlyphs[ n ] = NULL;
    mVolumeToRASSlice[ n ] = NULL;
    mSliceTransform[ n ] = NULL;
    mPlaneMappers[ n ] = NULL;
    mPlaneActors[ n ] = NULL;
    mGlyphEdgeActors[ n ] = NULL;
  }
    
}

vtkKWScubaLayer3DDTI::~vtkKWScubaLayer3DDTI () {
  
  mVolumeToRAS->Delete();
  mVolumeToRAS = NULL;

  mReducedVolume->Delete();
  mReducedVolume = NULL;

  for( int n=0; n<3; n++ ) {

    if ( mClips[ n ] )
      mClips[ n ]->Delete();

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
vtkKWScubaLayer3DDTI::SetDTIProperties ( ScubaCollectionPropertiesDTI* const iProperties ) {
  mDTIProperties = iProperties;
}


void
vtkKWScubaLayer3DDTI::Create () {

  // Bail if we don't have our source yet.
  if( NULL == mDTIProperties )
    throw runtime_error( "vtkKWScubaLayer3DDTI::Create: No source" );

  // Get some values from the MRI.
  float RASBounds[6];
  this->GetUnscaledRASBounds( RASBounds );

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

  // This rotates the volume to the proper orientation.
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
  
  // this is used for our fast rendering
  mReducedVolume = vtkImageShrink3D::New();
  mReducedVolume->SetInputConnection( mVolumeToRAS->GetOutputPort() );

  // set up the transform for reorienting the tensors
  double *vtr = mDTIProperties->GetFAVolumeSource()->GetVoxelToRASMatrix();
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

  // three times for three cut planes
  for( int n = 0; n<3; n++ ) {  

    mClips[ n ] = vtkImageClip::New();

    // we don't want to display everything, so let's not displaying anything at 
    // the moment, and leave it to UpdatePlanes() to set the plane to display
    mClips[ n ]->SetOutputWholeExtent( 0, 0, 0, 0, 0, 0 );
    mClips[ n ]->SetInputConnection( mReducedVolume->GetOutputPort() );

    mGlyphs[ n ] = vtkFDTensorGlyph::New();
    mGlyphs[ n ]->SetInputConnection( mClips[ n ]->GetOutputPort() );

    // this reorients each tensor because when we resample the volume, the vector
    // data isn't reoriented
    mGlyphs[ n ]->SetVoxelToMeasurementFrameTransform( voxelTransform );

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

      this->PipelineChanged();
    }

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
vtkKWScubaLayer3DDTI::GetUnscaledRASBounds ( float ioBounds[6] ) const {
    
  this->GetRASBounds( ioBounds );
  
  const float scales[] = {
    mDTIProperties->GetFAVolumeSource()->GetPixelSizeX(),
    mDTIProperties->GetFAVolumeSource()->GetPixelSizeY(),
    mDTIProperties->GetFAVolumeSource()->GetPixelSizeZ()
  };
    
  // unscale the bounds
  for( int nDim=0; nDim<3; nDim++ ) {
    
    // unscale the bottom
    int nBound = nDim * 2;    
    ioBounds[ nBound ] = ioBounds[ nBound ] / scales[ nDim ];    
    
    // unscale the top
    nBound++;
    ioBounds[ nBound ] = ioBounds[ nBound ] / scales[ nDim ];
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
  
  const int ras[ 3 ] = {
    static_cast< int >( mViewProperties->Get3DRASX() / mDTIProperties->GetFAVolumeSource()->GetPixelSizeX() ),
    static_cast< int >( mViewProperties->Get3DRASY() / mDTIProperties->GetFAVolumeSource()->GetPixelSizeY() ),
    static_cast< int >( mViewProperties->Get3DRASZ() / mDTIProperties->GetFAVolumeSource()->GetPixelSizeZ() )
  };

  float floatBounds[ 6 ];
  this->GetUnscaledRASBounds ( floatBounds );
  
  // cast them as ints so that we don't get warnings
  int bounds[ 6 ];
  for( int cBound=0; cBound<6; cBound++ ) {
    bounds[ cBound ] = static_cast< int >( floatBounds[ cBound ] );
  }

  if( mIsShowingAll ) {

    mClips[ 0 ]->SetOutputWholeExtent( 
      bounds[ 0 ], bounds[ 1 ], 
      bounds[ 2 ], bounds[ 3 ], 
      bounds[ 4 ], bounds[ 5 ] );

  } else {
  
    const int currentShrinkage = this->GetCurrentShrinkage();
      
    // only update the slice that changed
    if( bounds[ 0 ] != ras[ 0 ] || bounds[ 1 ] != ras[ 0 ] || hasDetailChanged ) {
      
      // the current plane that we're modifying
      const int dim = 0;
      
      // moves the tensor plane up a little, so that it's cut in half by the 3D
      // mri and well as moving it up a little whenever the RAS moves up, but
      // the tensor slice volume isn't necessarily at a new slice yet
      const float reposition = ( static_cast< int >( ras[ dim ] ) % static_cast< int >( currentShrinkage ) ) + 1;
      this->ResetSliceTransform( dim );        
      mSliceTransform[dim]->Translate( reposition, 0, 0 );

      // current slice position to extract from the tensor volume
      const int position = static_cast< int >( floor( (float)(ras[ dim ] / currentShrinkage) ) );
      mClips[ dim ]->SetOutputWholeExtent( 
        position, position, 
        bounds[ 2 ], bounds[ 3 ], 
        bounds[ 4 ], bounds[ 5 ] );
    }
  
    if( bounds[ 2 ] != ras[ 1 ] || bounds[ 3 ] != ras[ 1 ] || hasDetailChanged ) {

      const int dim = 1;
      
      const float reposition = ( static_cast< int >( ras[ dim ] ) % static_cast< int >( currentShrinkage ) ) + 1;
      this->ResetSliceTransform( dim );
      mSliceTransform[dim]->Translate( 0, reposition, 0 );

      const int position = static_cast< int >( floor( (float)(ras[ dim ] / currentShrinkage) ) );
      mClips[ dim ]->SetOutputWholeExtent( 
        bounds[ 0 ], bounds[ 1 ], 
        position, position, 
        bounds[ 4 ], bounds[ 5 ] );
    }
  
    if( bounds[ 4 ] != ras[ 2 ] || bounds[ 5 ] != ras[ 2 ] || hasDetailChanged ) {
      
      const int dim = 2;
      
      const float reposition = ( static_cast< int >( ras[ dim ] ) % static_cast< int >( currentShrinkage ) ) + 1;
      this->ResetSliceTransform( dim );
      mSliceTransform[dim]->Translate( 0, 0, reposition );

      const int position = static_cast< int >( floor( (float)(ras[ dim ] / currentShrinkage) ) );
      mClips[ dim ]->SetOutputWholeExtent( 
        bounds[ 0 ], bounds[ 1 ],
        bounds[ 2 ], bounds[ 3 ],
        position, position );
        
    }

  }
  
  this->PipelineChanged();
  
}

void 
vtkKWScubaLayer3DDTI::ResetSliceTransform( const int iDim ) {
 
  mSliceTransform[ iDim ]->Identity();
  
  mSliceTransform[ iDim ]->Scale( 
    mDTIProperties->GetFAVolumeSource()->GetPixelSizeX(),
    mDTIProperties->GetFAVolumeSource()->GetPixelSizeY(),
    mDTIProperties->GetFAVolumeSource()->GetPixelSizeZ() );  

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
  
  if( mReducedVolume ) {
    
    // figure out the level of detail
    int shrinkage = this->GetCurrentShrinkage();
        
    for( int n = 0; n<3; n++ ) {          
      mReducedVolume->SetShrinkFactors( shrinkage, shrinkage, shrinkage );
    }
    mReducedVolume->AveragingOn();

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
    shrinkage = 6;
    break;
  case ScubaCollectionPropertiesDTI::Less:
    shrinkage = 3;
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
