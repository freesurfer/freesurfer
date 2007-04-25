/**
 * @file  vtkKWScubaLayer3DPath.cxx
 * @brief Displays 3D path
 *
 */
/*
 * Original Author: Dennis Jen
 * CVS Revision Info:
 *    $Author: dsjen $
 *    $Date: 2007/04/25 19:24:48 $
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

#include "vtkKWRadioButtonSet.h"
#include "vtkKWRadioButtonSetWithLabel.h"
#include "vtkKWRadioButton.h"
#include "vtkKWScaleWithEntry.h"
#include "vtkKWCheckButton.h"

using namespace std;

vtkStandardNewMacro( vtkKWScubaLayer3DPath );
vtkCxxRevisionMacro( vtkKWScubaLayer3DPath, "$Revision: 1.4 $" );

vtkKWScubaLayer3DPath::vtkKWScubaLayer3DPath () :
  mPathProperties( NULL ),
  mPathActor( NULL ),
  mStartPointActor( NULL ),
  mEndPointActor( NULL ),
  mPathPointsCollection( NULL ),
  mTubeFilter( NULL ),
  mTubeMapper( NULL ),
  mTubeActor( NULL ),
  mScaleTubeRadius( NULL ),
  mTriangleStripper( NULL ),
  mNormals( NULL ),
  mMapper( NULL ),
  mPathMode( TUBE_MODE ),
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
    
  if( mTriangleStripper ) {
    mTriangleStripper->Delete();
  }
  
  if( mNormals ) {
    mNormals->Delete();
  }
  
  if( mMapper ) {
    mMapper->Delete();
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
vtkKWScubaLayer3DPath::SetPathProperties ( ScubaCollectionPropertiesPath* const iProperties ) {
  mPathProperties = iProperties;
}


void
vtkKWScubaLayer3DPath::Create () {

  // Bail if we don't have our source yet.
  if( NULL == mPathProperties )
    throw runtime_error( "vtkKWScubaLayer3DPath::Create: No source" );

  // this creates the tube representation of the path    
  this->CreateTube();
    
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
  radBtnSetPathMode->GetWidget()->PackHorizontallyOn();
  radBtnSetPathMode->SetLabelText( "Mode: " );
  
  char sTclCmd[1024];

  vtkKWRadioButton* radBtnTubeMode = radBtnSetPathMode->GetWidget()->AddWidget( 0 );
  radBtnTubeMode->SetText( "Tube" );
  sprintf( sTclCmd, "SetPathMode %d", TUBE_MODE );
  radBtnTubeMode->SetCommand( this, sTclCmd );
  if ( mPathMode == TUBE_MODE )
    radBtnTubeMode->SelectedStateOn();

  vtkKWRadioButton* radBtnThresholdMode = radBtnSetPathMode->GetWidget()->AddWidget( 1 );
  radBtnThresholdMode->SetText( "Threshold" );
  sprintf( sTclCmd, "SetPathMode %d", THRESHOLD_MODE );
  radBtnThresholdMode->SetCommand( this, sTclCmd );
  if ( mPathMode == THRESHOLD_MODE )
    radBtnThresholdMode->SelectedStateOn();

  radBtnTubeMode->Delete();
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
  vtkPolyData *polyData = vtkPolyData::New();
  polyData->SetPoints( inputPoints );
  polyData->SetLines( lines );

  // set up the scalar data that was sampled by the pathway
  vtkFloatArray *sampleData = vtkFloatArray::New();
  sampleData->SetName( "SampleData" );
  
  // error check to make sure that number of samples is the same as the number
  // of points
  if( nPoints == mPathProperties->GetNumberOfSamples() ) {
    // get the scalar data
    for( int i=0; i<mPathProperties->GetNumberOfSamples(); i++ ) {
      const double value = mPathProperties->GetPointSampleValue( i );
      sampleData->InsertNextValue( value );
    }
    polyData->GetPointData()->SetScalars( sampleData );
  }
  
  inputPoints->Delete();
  lines->Delete();
  sampleData->Delete();
  
  // Add thickness to the resulting line.
  mTubeFilter = vtkTubeFilter::New();
  mTubeFilter->SetNumberOfSides( 5 );
  mTubeFilter->SetInput( polyData );
  mTubeFilter->SetRadius( mTubeRadius );
  
  if( mbIsTubeRadiusScaled ) {
    mTubeFilter->SetVaryRadiusToVaryRadiusByScalar();
  } else {
    mTubeFilter->SetVaryRadiusToVaryRadiusOff();
  }
      
  polyData->Delete();
  
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
    
    if( mPathMode == TUBE_MODE ) {
      
      // turn the threshold visualization off and the tube on      
      if( mPathActor ) {
        mPathActor->SetVisibility( false );
      }
      
      if( mTubeActor ) {
        mTubeActor->SetVisibility( true );
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
