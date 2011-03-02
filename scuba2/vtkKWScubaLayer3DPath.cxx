/**
 * @file  vtkKWScubaLayer3DPath.cxx
 * @brief Displays 3D path
 *
 */
/*
 * Original Author: Dennis Jen
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:40 $
 *    $Revision: 1.12 $
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


#include "vtkKWScubaLayer3DPath.h"

#include <string>
#include <stdexcept>

#include "ScubaCollectionProperties.h"
#include "ScubaCollectionPropertiesPath.h"
#include "ScubaInfoItem.h"

#include "vtkFSVolumeSource.h"

#include "vtkActor.h"
#include "vtkActorCollection.h"
#include "vtkObjectFactory.h"

#include "vtkPolyDataMapper.h"
#include "vtkProperty.h"
#include "vtkPolyDataNormals.h"
#include "vtkStripper.h"
#include "vtkSimplePointsReader.h"
#include "vtkSphereSource.h"
#include "vtkPoints.h"
#include "vtkCellArray.h"
#include "vtkTubeFilter.h"
#include "vtkFloatArray.h"
#include "vtkPointData.h"
#include "vtkTransform.h"

#include "vtkKWRadioButtonSet.h"
#include "vtkKWRadioButtonSetWithLabel.h"
#include "vtkKWRadioButton.h"
#include "vtkKWScaleWithEntry.h"
#include "vtkKWCheckButton.h"

using namespace std;

vtkStandardNewMacro( vtkKWScubaLayer3DPath );
vtkCxxRevisionMacro( vtkKWScubaLayer3DPath, "$Revision: 1.12 $" );

vtkKWScubaLayer3DPath::vtkKWScubaLayer3DPath () :
  mPathProperties( NULL ),
  mPathActor( NULL ),
  mStartPointActor( NULL ),
  mEndPointActor( NULL ),
  mPathPointsCollection( NULL ),
  mInitialPathsActorCollection( NULL ),
  mTubePolyData( NULL ),
  mSampleData( NULL ),
  mProbabilityData( NULL ),
  mTubeFilter( NULL ),
  mTubeMapper( NULL ),
  mTubeActor( NULL ),
  mScaleTubeRadius( NULL ),
  mTriangleStripper( NULL ),
  mNormals( NULL ),
  mMapper( NULL ),
  mPathMode( PROBABILITY_TUBE_MODE ),
  mTubeRadius( 1.0 ),
  mbIsTubeRadiusScaled( true ),
  mbIsTubeRadiusColored( true ) {
};

vtkKWScubaLayer3DPath::~vtkKWScubaLayer3DPath () {
    
  if( mPathActor ) {
    mPathActor->Delete();
  }
  
  if( mStartPointActor ) {
    mStartPointActor->Delete();
  }
  
  if( mEndPointActor ) {
    mEndPointActor->Delete();
  }

  if( mPathPointsCollection ) {
    mPathPointsCollection->Delete();
  }
  
  if( mInitialPathsActorCollection ) {
    mInitialPathsActorCollection->Delete();
  }
    
  if( mTriangleStripper ) {
    mTriangleStripper->Delete();
  }
  
  if( mNormals ) {
    mNormals->Delete();
  }
  
  if( mMapper ) {
    mMapper->Delete();
  }

  if( mTubePolyData ) {
    mTubePolyData->Delete();
  }
  
  if( mSampleData ) {
    mSampleData->Delete();
  }
  
  if( mProbabilityData ) {
    mProbabilityData->Delete();
  }
  
  if( mTubeFilter ) {
    mTubeFilter->Delete();
  }
  
  if( mTubeMapper ) {
    mTubeMapper->Delete();
  }
  
  if( mTubeActor ) {
    mTubeActor->Delete();
  }
  
  if( mScaleTubeRadius ) {
    mScaleTubeRadius->Delete();
  }
  
}

void
vtkKWScubaLayer3DPath::SetPathProperties ( ScubaCollectionPropertiesPath const* iProperties ) {
  mPathProperties = iProperties;
}


void
vtkKWScubaLayer3DPath::Create () {

  // Bail if we don't have our source yet.
  if( NULL == mPathProperties )
    throw runtime_error( "vtkKWScubaLayer3DPath::Create: No source" );

  // this creates the tube representation of the path    
  this->CreateTube();
  
  // create the multiple initial paths
  this->CreateInitialTubes();
    
  mNormals = vtkPolyDataNormals::New();
  mNormals->SetInput( mPathProperties->GetMesh() );

  // create triangle strips for faster rendering
  mTriangleStripper = vtkStripper::New();
  mTriangleStripper->SetInputConnection( mNormals->GetOutputPort() );
  
  // set up the mapper and don't color it by scalars or it'll be all red
  mMapper = vtkPolyDataMapper::New();
  mMapper->SetInputConnection( mTriangleStripper->GetOutputPort() );
  mMapper->ScalarVisibilityOff();
  
  // set up the actor for the pathway
  mPathActor = vtkActor::New();
  mPathActor->SetMapper( mMapper );
  mPathActor->GetProperty()->SetAmbient( 0.5 );
  mPathActor->GetProperty()->SetDiffuse( 0.2 );
  mPathActor->GetProperty()->SetSpecular( 0.2 );
    
  this->UpdatePathColor();
    
  this->UpdateOpacity();
  
  this->UpdatePathMode();
  
  this->AddProp( mPathActor );

//// TODO: 
  // display the start point
//  const double radius = 3.0;
//  const int startPoint = 0;
//  mStartPointActor = vtkActor::New();
//  this->CreatePoint( mStartPointActor, startPoint, radius );
//  
//  // display the end point
//  const int endPoint = 
//    mPathProperties->GetPathPointsSource()->GetOutput()->GetNumberOfPoints() - 1;
//  mEndPointActor = vtkActor::New();
//  this->CreatePoint( mEndPointActor, endPoint, radius );
  
//  mPathPointsCollection = vtkActorCollection::New();
//  const int nPoints = mPathProperties->GetPathPointsSource()->GetOutput()->GetNumberOfPoints();
//  for( int i=0; i<nPoints; i++ ) {
//
//    // the actor must be created here, rather than in the CreatePoint method.
//    // I think it has something to do with memory management and the actor being
//    // deleted otherwise.
//    vtkActor *actor = vtkActor::New();
//    const float radius = 1.0;
//    this->CreatePoint( actor, i, radius );
//            
//    mPathPointsCollection->AddItem( actor );
//    
//  }
  
}

void
vtkKWScubaLayer3DPath::DoListenToMessage ( string isMessage, void* iData ) {

  if( isMessage == "OpacityChanged" ) {
    this->UpdateOpacity();
  } else if (isMessage == "PathColorChanged" ) {
    this->UpdatePathColor();
  }
  
}

void 
vtkKWScubaLayer3DPath::SetPathMode( int iMode ) {
  if( iMode != mPathMode ) {
    mPathMode = iMode;
    this->UpdatePathMode();
  }
}

void
vtkKWScubaLayer3DPath::SetTubeRadius( float iRadius ) {
  
  if( mTubeRadius != iRadius ) {
    
    mTubeRadius = iRadius;
  
    if( mTubeFilter ) {    
      mTubeFilter->SetRadius( mTubeRadius );    
      this->PipelineChanged();
    }
    
  }
  
}

void 
vtkKWScubaLayer3DPath::SetScaleTubeRadius( int ibIsScaling ) {

  if( ibIsScaling != mbIsTubeRadiusScaled ) {
    mbIsTubeRadiusScaled = ibIsScaling;
    
    if( mTubeFilter ) {
    
      if( mbIsTubeRadiusScaled ) {
        mTubeFilter->SetVaryRadiusToVaryRadiusByScalar();
      } else {
        mTubeFilter->SetVaryRadiusToVaryRadiusOff();
      }
      
      this->PipelineChanged();
    
    }
    
  }
  
}

void 
vtkKWScubaLayer3DPath::SetColoredTube( int ibIsColored ) {
  
  if( ibIsColored != mbIsTubeRadiusColored ) {
    mbIsTubeRadiusColored = ibIsColored;

    if( mTubeMapper ) {
    
      if( mbIsTubeRadiusColored ) {
        mTubeMapper->ScalarVisibilityOn();
      } else {
        mTubeMapper->ScalarVisibilityOff();
      }
      
      this->PipelineChanged();
    
    }
    
  }
  
}


void 
vtkKWScubaLayer3DPath::AddControls ( vtkKWWidget* iPanel ) {
  
  // Path Representation
  vtkKWRadioButtonSetWithLabel* radBtnSetPathMode = vtkKWRadioButtonSetWithLabel::New();
  radBtnSetPathMode->SetParent( iPanel );
  radBtnSetPathMode->Create();
  
  // in other words, pack it vertically
  radBtnSetPathMode->GetWidget()->PackHorizontallyOff();
  radBtnSetPathMode->SetLabelText( "Mode: " );
  
  char sTclCmd[1024];

  vtkKWRadioButton* radBtnProbabilityTubeMode = radBtnSetPathMode->GetWidget()->AddWidget( 0 );
  radBtnProbabilityTubeMode->SetText( "Probability Tube" );
  sprintf( sTclCmd, "SetPathMode %d", PROBABILITY_TUBE_MODE );
  radBtnProbabilityTubeMode->SetCommand( this, sTclCmd );
  if ( mPathMode == PROBABILITY_TUBE_MODE )
    radBtnProbabilityTubeMode->SelectedStateOn();

  vtkKWRadioButton* radBtnSampledTubeMode = radBtnSetPathMode->GetWidget()->AddWidget( 1 );
  radBtnSampledTubeMode->SetText( "Sampled Tube" );
  sprintf( sTclCmd, "SetPathMode %d", SAMPLED_TUBE_MODE );
  radBtnSampledTubeMode->SetCommand( this, sTclCmd );
  if ( mPathMode == SAMPLED_TUBE_MODE )
    radBtnSampledTubeMode->SelectedStateOn();

  vtkKWRadioButton* radBtnThresholdMode = radBtnSetPathMode->GetWidget()->AddWidget( 2 );
  radBtnThresholdMode->SetText( "Threshold" );
  sprintf( sTclCmd, "SetPathMode %d", THRESHOLD_MODE );
  radBtnThresholdMode->SetCommand( this, sTclCmd );
  if ( mPathMode == THRESHOLD_MODE )
    radBtnThresholdMode->SelectedStateOn();
    
  radBtnProbabilityTubeMode->Delete();
  radBtnSampledTubeMode->Delete();
  radBtnThresholdMode->Delete();
  
  // create a widget for adjusting the radius of the tube
  mScaleTubeRadius = vtkKWScaleWithEntry::New();
  mScaleTubeRadius->SetParent( iPanel );
  mScaleTubeRadius->SetOrientationToHorizontal();
  mScaleTubeRadius->Create();
  mScaleTubeRadius->SetLabelText( "Tube Radius: " );
  mScaleTubeRadius->SetRange( 0, 5 );
  mScaleTubeRadius->SetResolution( 0.1 );
  mScaleTubeRadius->SetEntryWidth( 3 );
  mScaleTubeRadius->SetCommand( this, "SetTubeRadius" );
  mScaleTubeRadius->SetValue( mTubeRadius );
  
  // check box for scaling radius of tube by sampled scalar
  mChkBtnScaleTubeRadius = vtkKWCheckButton::New();
  mChkBtnScaleTubeRadius->SetParent( iPanel );
  mChkBtnScaleTubeRadius->Create();
  mChkBtnScaleTubeRadius->SetAnchorToWest();
  mChkBtnScaleTubeRadius->SetText( "Scale Tube Radius" );
  mChkBtnScaleTubeRadius->SetCommand( this, "SetScaleTubeRadius" );
  if ( mbIsTubeRadiusScaled )
    mChkBtnScaleTubeRadius->SelectedStateOn();
    
  // check box for coloring tube by sampled scalar
  mChkBtnColorTubeRadius = vtkKWCheckButton::New();
  mChkBtnColorTubeRadius->SetParent( iPanel );
  mChkBtnColorTubeRadius->Create();
  mChkBtnColorTubeRadius->SetAnchorToWest();
  mChkBtnColorTubeRadius->SetText( "Color Tube Radius" );
  mChkBtnColorTubeRadius->SetCommand( this, "SetColoredTube" );
  if ( mbIsTubeRadiusColored )
    mChkBtnColorTubeRadius->SelectedStateOn();
  
  this->Script( "pack %s %s %s %s -side top -fill x -anchor nw",
                radBtnSetPathMode->GetWidgetName(),
                mScaleTubeRadius->GetWidgetName(),
                mChkBtnScaleTubeRadius->GetWidgetName(),
                mChkBtnColorTubeRadius->GetWidgetName() );
                
}

void 
vtkKWScubaLayer3DPath::RemoveControls () {
  
  if( mScaleTubeRadius ) {
    mScaleTubeRadius->Delete();
    mScaleTubeRadius = NULL;
  }
  
  if( mChkBtnScaleTubeRadius ) {
    mChkBtnScaleTubeRadius->Delete();
    mChkBtnScaleTubeRadius = NULL;
  }
  
  if( mChkBtnColorTubeRadius ) {
    mChkBtnColorTubeRadius->Delete();
    mChkBtnColorTubeRadius = NULL;
  }
  
}

void 
vtkKWScubaLayer3DPath::CreatePoint ( vtkActor* pointActor, const int iPointIndex, const double iRadius ) {
  
  vtkFSVolumeSource* source = mPathProperties->GetPathVolumeSource();
  
  float RASX, RASY, RASZ;
  double *point = 
    mPathProperties->GetPathPointsSource()->GetOutput()->GetPoint( iPointIndex );
  
  source->ConvertIndexToRAS( point[ 0 ], point[ 1 ], point[ 2 ], 
    RASX, RASY, RASZ );
    
  // the center of the bounding volume is actually in the lower corner, so we
  // need to move our transformed points diagonally lower by the amount of the
  // world's center
  float worldCenter[3];
  worldCenter[0] = mPathProperties->GetPathVolumeSource()->GetRASCenterX();
  RASX -= worldCenter[0];

  worldCenter[1] = mPathProperties->GetPathVolumeSource()->GetRASCenterY();
  RASY -= worldCenter[1];

  worldCenter[2] = mPathProperties->GetPathVolumeSource()->GetRASCenterZ();
  RASZ -= worldCenter[2];
  
  vtkSphereSource* sphere = vtkSphereSource::New();
  sphere = vtkSphereSource::New();
  sphere->SetCenter( RASX, RASY, RASZ );
  sphere->SetRadius( iRadius );

  vtkPolyDataMapper* mapper = vtkPolyDataMapper::New();
  mapper->SetInputConnection( sphere->GetOutputPort() );

  pointActor->SetMapper( mapper );
  
  double r, g, b;
  mPathProperties->GetPointColor( iPointIndex, r, g, b );
  pointActor->GetProperty()->SetColor( r, g, b );

  this->AddProp( pointActor );

}

void
vtkKWScubaLayer3DPath::CreateTube() {
  
  const int nPoints = mPathProperties->GetPathPointsSource()->GetOutput()->GetNumberOfPoints();
  
  vtkPoints *inputPoints = vtkPoints::New();
  vtkCellArray *lines = vtkCellArray::New();
  lines->InsertNextCell( nPoints );
  
  vtkFSVolumeSource* source = mPathProperties->GetPathVolumeSource();

  float worldCenter[3];
  worldCenter[0] = source->GetRASCenterX();  
  worldCenter[1] = source->GetRASCenterY();
  worldCenter[2] = source->GetRASCenterZ();
  
  // insert the tube data
  for( int i=0; i<nPoints; i++ ) {

    double *point = 
      mPathProperties->GetPathPointsSource()->GetOutput()->GetPoint( i );

    float RASX, RASY, RASZ;
    source->ConvertIndexToRAS( point[ 0 ], point[ 1 ], point[ 2 ], 
      RASX, RASY, RASZ );

    // the center of the bounding volume is actually in the lower corner, so we
    // need to move our transformed points diagonally lower by the amount of the
    // world's center
    RASX -= worldCenter[0];
    RASY -= worldCenter[1];
    RASZ -= worldCenter[2];
    
    inputPoints->InsertPoint( i, RASX, RASY, RASZ );
    lines->InsertCellPoint( i );    
  }
  
  // add the points, lines, and sample data    
  mTubePolyData = vtkPolyData::New();
  mTubePolyData->SetPoints( inputPoints );
  mTubePolyData->SetLines( lines );

  // set up the probablity data
  mProbabilityData = vtkFloatArray::New();
  mProbabilityData->SetName( "ProbablityData" );
  
  // error check to make sure that number of samples is the same as the number
  // of points
  if( nPoints == mPathProperties->GetNumberOfProbabilities() ) {
    // get the scalar data
    for( int i=0; i<mPathProperties->GetNumberOfProbabilities(); i++ ) {
      const double value = mPathProperties->GetPointProbabilityValue( i );
      mProbabilityData->InsertNextValue( value );
    }
  }

  // set up the scalar data that was sampled by the pathway
  mSampleData = vtkFloatArray::New();
  mSampleData->SetName( "SampleData" );
  
  // error check to make sure that number of samples is the same as the number
  // of points
  if( nPoints == mPathProperties->GetNumberOfSamples() ) {
    // get the scalar data
    for( int i=0; i<mPathProperties->GetNumberOfSamples(); i++ ) {
      const double value = mPathProperties->GetPointSampleValue( i );
      mSampleData->InsertNextValue( value );
    }
  }
  
  // set the right data as the scalars initially
  if( mPathMode == PROBABILITY_TUBE_MODE ) {
     mTubePolyData->GetPointData()->SetScalars( mProbabilityData );        
  } else if( mPathMode == SAMPLED_TUBE_MODE ) {
     mTubePolyData->GetPointData()->SetScalars( mSampleData );        
  }
  
  
  inputPoints->Delete();
  lines->Delete();
  
  // Add thickness to the resulting line.
  mTubeFilter = vtkTubeFilter::New();
  mTubeFilter->SetNumberOfSides( 5 );
  mTubeFilter->SetInput( mTubePolyData );
  mTubeFilter->SetRadius( mTubeRadius );
  
  if( mbIsTubeRadiusScaled ) {
    mTubeFilter->SetVaryRadiusToVaryRadiusByScalar();
  } else {
    mTubeFilter->SetVaryRadiusToVaryRadiusOff();
  }
    
  mTubeMapper = vtkPolyDataMapper::New();
  mTubeMapper->SetInputConnection( mTubeFilter->GetOutputPort() );
  
  if( mbIsTubeRadiusColored ) {
    mTubeMapper->ScalarVisibilityOn();
  } else {
    mTubeMapper->ScalarVisibilityOff();
  }

  mTubeActor = vtkActor::New();
  mTubeActor->SetMapper( mTubeMapper );
  
  this->AddProp( mTubeActor );
}

void 
vtkKWScubaLayer3DPath::CreateInitialTubes() {
  
  // get the initial paths
  const std::vector< std::vector< double* >* > * initialPaths = mPathProperties->GetInitialPaths();  
  
  // only create the initial paths if they're available...
  if( initialPaths != NULL && !initialPaths->empty() ) {
    
    vtkFSVolumeSource* source = mPathProperties->GetPathVolumeSource();
    
    float worldCenter[3];
    worldCenter[0] = source->GetRASCenterX();  
    worldCenter[1] = source->GetRASCenterY();
    worldCenter[2] = source->GetRASCenterZ();
    
    // this is a long for loop to iterate over the initial paths
    for( std::vector< std::vector< double* >* >::const_iterator itPaths = initialPaths->begin(); 
      itPaths != initialPaths->end(); itPaths++ ) {
        
      std::vector< double* >* path = *itPaths;

      vtkPoints *inputPoints = vtkPoints::New();
      vtkCellArray *lines = vtkCellArray::New();
      lines->InsertNextCell( path->size() );
      
      int nPoint = 0;
      // another long for loop to iterate over the points in a path
      for( std::vector< double* >::const_iterator itPoints = path->begin(); 
        itPoints != path->end(); itPoints++ ) {
          
        double *point = ( *itPoints );
        
        // convert each point
        float RASX, RASY, RASZ;
        source->ConvertIndexToRAS( point[ 0 ], point[ 1 ], point[ 2 ], 
          RASX, RASY, RASZ );
    
        // the center of the bounding volume is actually in the lower corner, so we
        // need to move our transformed points diagonally lower by the amount of the
        // world's center
        RASX -= worldCenter[0];
        RASY -= worldCenter[1];
        RASZ -= worldCenter[2];
        
        inputPoints->InsertPoint( nPoint, RASX, RASY, RASZ );
        lines->InsertCellPoint( nPoint );
        
        nPoint++;
        
      }
      
      // add the points and lines
      vtkPolyData *polyData = vtkPolyData::New();
      polyData->SetPoints( inputPoints );
      polyData->SetLines( lines );      
      
      inputPoints->Delete();
      lines->Delete();
      
      // Add thickness to the resulting line.
      vtkTubeFilter *tubeFilter = vtkTubeFilter::New();
      tubeFilter->SetNumberOfSides( 5 );
      tubeFilter->SetInput( polyData );
      tubeFilter->SetRadius( 0.5 );
      tubeFilter->SetVaryRadiusToVaryRadiusOff();
          
      polyData->Delete();
      
      vtkPolyDataMapper *tubeMapper = vtkPolyDataMapper::New();
      tubeMapper->SetInputConnection( tubeFilter->GetOutputPort() );      
      tubeMapper->ScalarVisibilityOff();
    
      vtkActor *tubeActor = vtkActor::New();
      tubeActor->SetMapper( tubeMapper );
      
      // pseudo randomizers colors
      const double r = ( rand() % 100 + 1 ) / 100.0;
      const double g = ( rand() % 100 + 1 ) / 100.0;
      const double b = 0.8;
      tubeActor->GetProperty()->SetColor( r, g, b );
      
      // TODO: add this to the collection too...
      this->AddProp( tubeActor );
      
    }
    
  }
    
}

void 
vtkKWScubaLayer3DPath::UpdateOpacity() {

  if( NULL != mPathProperties ) {
  
    if( mPathActor ) {
      mPathActor->GetProperty()->SetOpacity( mProperties->GetOpacity() );
    }

    if( mTubeActor ) {
      mTubeActor->GetProperty()->SetOpacity( mProperties->GetOpacity() );
    }

    this->PipelineChanged();
    
  }
  
}

void 
vtkKWScubaLayer3DPath::UpdatePathColor() {

  if( NULL != mPathProperties ) {
    
    double r, g, b;
    mPathProperties->GetPathColor( r, g, b );

    if( mPathActor ) {
      mPathActor->GetProperty()->SetColor( r, g, b );
    }
    
    if( mTubeActor ) {
      mTubeActor->GetProperty()->SetColor( r, g, b );
    }
    
  }
  
}

void 
vtkKWScubaLayer3DPath::UpdatePathMode() {

  if( NULL != mPathProperties ) {
    
    if( mPathMode == PROBABILITY_TUBE_MODE || mPathMode == SAMPLED_TUBE_MODE ) {
      
      // turn the threshold visualization off and the tube on      
      if( mPathActor ) {
        mPathActor->SetVisibility( false );
      }
      
      if( mTubeActor ) {        
        mTubeActor->SetVisibility( true );
      }
      
      if( mPathMode == PROBABILITY_TUBE_MODE ) {
         mTubePolyData->GetPointData()->SetScalars( mProbabilityData );        
      } else if( mPathMode == SAMPLED_TUBE_MODE ) {
         mTubePolyData->GetPointData()->SetScalars( mSampleData );        
      }
      
    } else if (  mPathMode == THRESHOLD_MODE ) {

      // turn the tube off and the threshold on
      if( mPathActor ) {
        mPathActor->SetVisibility( true );
      }
      
      if( mTubeActor ) {
        mTubeActor->SetVisibility( false );
      }
      
    }
    
    this->PipelineChanged();
    
  }
  
}
