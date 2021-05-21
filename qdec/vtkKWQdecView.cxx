/**
 * @brief The VTK pipeline and Renderer
 *
 * Builds the VTK pipeline to display the surface and scalars. Also
 * manages the callback for clicking on a vertex and passing the point
 * along to the Fsgdf plotting code.
 */
/*
 * Original Author: Kevin Teich
 *
 * Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
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

#include <stdexcept>
#include <sstream>
#include <assert.h>
#include "vtkKWQdecView.h"

#include "QdecEvents.h"
#include "QdecUtilities.h"
#include "QdecVertexAnnotationLookup.h"
#include "vtkAbstractPicker.h"
#include "vtkActor.h"
#include "vtkCamera.h"
#include "vtkCellArray.h"
#include "vtkCellPicker.h"
#include "vtkColorTransferFunction.h"
#include "vtkCommand.h"
#include "vtkCornerAnnotation.h"
#include "vtkCursor3D.h"
#include "vtkDataArray.h"
#include "vtkFloatArray.h"
#include "vtkFSSurfaceSource.h"
#include "vtkFSSurfaceScalarsReader.h"
#include "vtkKWApplication.h"
#include "vtkKWProgressDialog.h"
#include "vtkInflatePolyData.h"
#include "vtkLODProp3D.h"
#include "vtkMatrix4x4.h"
#include "vtkObjectFactory.h"
#include "vtkOutlineFilter.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkPolyDataMapper.h"
#include "vtkProp3DCollection.h"
#include "vtkProperty.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkRenderer.h"
#include "vtkScalarBarActor.h"
#include "vtkTextProperty.h"
#include "vtkTransform.h"
#include "vtkTransformPolyDataFilter.h"

using namespace std;

vtkStandardNewMacro( vtkKWQdecView );

// these control the amount and speed of rotation
// with AnimateSteps=1, it doesnt animate, and its instaneous
// old values shown for the slow rotation animation when Rotate buttons pushed
double const vtkKWQdecView::kAnimationDegrees = 90.0;
double const vtkKWQdecView::kcAnimateSteps = 1;  // was 10.0
double const vtkKWQdecView::kAnimateTimeInSeconds = 0.1;  // was 1.0

vtkKWQdecView::vtkKWQdecView () :
  mfnSurface( "" ),
  mVertexPlot( NULL ),
  mnCursorVertexIndex( -1 ),
  mOverlayOpacity( 0.3 ),
  mAnnotationLookup( NULL ),
  mbInSelection( false ),
  mnFirstVertexInPath( -1 ),
  mnLastVertexInPath( -1 ) {
  for( int n = 0; n < 3; n++ ) {
    mDefaultPosition[n] = 0;
    mDefaultFocalPoint[n] = 0;
    mDefaultViewUp[n] = 0;
  }
  mDefaultZoom = 1.2; // make the default just a little bigger than normal
  
  mSurfaceSelectionPoints = vtkPoints::New();
  mSurfaceSelectionLines = vtkCellArray::New();
  
}

vtkKWQdecView::~vtkKWQdecView () {
  
  this->GetRenderer()->RemoveAllViewProps();
}

void
vtkKWQdecView::CreateWidget () {

  this->Superclass::CreateWidget();

  // Set the picker.
  vtkSmartPointer<vtkCellPicker> picker =
    vtkSmartPointer<vtkCellPicker>::New();
  this->GetRenderWindowInteractor()->SetPicker( picker );

  // Create our InteractorStyle and give it to the RWI.
  vtkSmartPointer<ViewInteractor> interactor =
    vtkSmartPointer<ViewInteractor>::New();
  interactor->SetView( this );
  this->GetRenderWindowInteractor()->SetInteractorStyle( interactor );
}

vtkKWQdecView::ViewInteractor*
vtkKWQdecView::ViewInteractor::New () {
  return new vtkKWQdecView::ViewInteractor();
}

vtkKWQdecView::ViewInteractor::ViewInteractor () :
  mView( NULL ), mnButtonDown( 0 ) {
}

void
vtkKWQdecView::ViewInteractor::SetView( vtkKWQdecView* iView ) {
  mView = iView;
}

void
vtkKWQdecView::ViewInteractor::OnLeftButtonDown () {

  assert( mView.GetPointer() );

  vtkInteractorStyleTrackballCamera::OnLeftButtonDown();

  mnButtonDown = 1;

  // whether button-1 or ctrl-button-1 selects the vertex is configurable
  // in the .Qdecrc file.  if GDF_BUTTON=BUTTON_1 is in the .Qdecrc file,
  // then button-1 selects the vertex, else ctrl-button-1 does that
  const char* sGdfButton = 
    QdecUtilities::GetQdecrcResourceString( "GDF_BUTTON" );

  if( this->GetInteractor()->GetShiftKey() ) {
    
    // If we're drawing paths, mark this as the first and last
    // vertex in the path.
    int nVertex = this->GetVertexAtPicker();
    if( -1 == nVertex ) // We may not be able to find a point. It happens.
      return;

    mView->mbInSelection = true;
    mView->mnFirstVertexInPath = nVertex;
    mView->mnLastVertexInPath = nVertex;

    // Clear the path.
    mView->ResetPath();

  } else if( this->GetInteractor()->GetControlKey() ) {

    if (sGdfButton && (0 == strcmp(sGdfButton,"BUTTON_1"))) {
      // just let it rotate camera
      // (need this else if so that it doesnt select vertex)
    } else {
      // Select vertex
      mView->SelectSurfaceVertex( this->GetVertexAtPicker() ); 
    }
  } else { 
    if (sGdfButton && (0 == strcmp(sGdfButton,"BUTTON_1"))) {
      // Select vertex
      mView->SelectSurfaceVertex( this->GetVertexAtPicker() ); 
    } else {
      // just let it rotate camera      
    }
  }
}

void
vtkKWQdecView::ViewInteractor::OnLeftButtonUp () {

  assert( mView.GetPointer() );

  vtkInteractorStyleTrackballCamera::OnLeftButtonUp();

  mnButtonDown = 0;

  // Drawing paths.
  if( mView->mbInSelection ) {

    // Join the last vertex to the first vertex.
    vector<int> lVertices;
    try {
      mView->mCurrentSurfaceSource->
        FindPath( mView->mnLastVertexInPath, mView->mnFirstVertexInPath,
                  lVertices );
    }
    catch( exception& e ) {
      return;
    }

    // Add all the vertices to the path.
    while( !lVertices.empty() ) {
      
      // Get a vertex number.
      int nVertexToAdd = lVertices.back();
      lVertices.pop_back();
      
      // Add it to the path.
      mView->AddVertexToPath( nVertexToAdd );
    }
    
    // Rebuild this segment and draw it.
    mView->RebuildPathLine ();
    mView->Render();
    
    // Clear the first and last vertex.
    mView->mnFirstVertexInPath = mView->mnLastVertexInPath = -1;
    mView->mbInSelection = false;
  }

}
void
vtkKWQdecView::ViewInteractor::OnMouseMove () {

  assert( mView.GetPointer() );

  // Drawing paths.
  if( mView->mbInSelection ) {

    // Find the surface vertex at this point.
    int nVertex = this->GetVertexAtPicker();
    if( -1 == nVertex ) // We may not be able to find a point. It happens.
      return;

    if( nVertex != mView->mnLastVertexInPath ) {
      
      // Make a path from the last clicked vertex to this one.
      vector<int> lVertices;
      try {
        mView->mCurrentSurfaceSource->
          FindPath( mView->mnLastVertexInPath, nVertex, lVertices );
      }
      catch( exception& e ) {
        return;
      }
      
      // Add all the vertices to the list of lines.
      while( !lVertices.empty() ) {
	
        // Get a vertex number.
        int nVertexToAdd = lVertices.back();
        lVertices.pop_back();
	    
        // Add it to the path.
        mView->AddVertexToPath( nVertexToAdd );
      }

      // Rebuild this segment and draw it.
      mView->RebuildPathLine ();
      mView->Render();
    }

    // Save this as the last point.
    mView->mnLastVertexInPath = nVertex;
    
    // End of drawing paths.
  } else {
  
    // Note we don't do this if we drawing paths, so we don't move the
    // camera.
    vtkInteractorStyleTrackballCamera::OnMouseMove();
  }
}

int
vtkKWQdecView::ViewInteractor::GetVertexAtPicker () {

  assert( mView.GetPointer() );

  // Try to get our picker.
  vtkCellPicker* picker =
    vtkCellPicker::SafeDownCast( this->GetInteractor()->GetPicker() );
  assert( picker );
  
  // Use the picker to pick and find the world coords of the point
  // picked.
  double worldCoords[3] = { 0,0,0 };
  picker->Pick( this->GetInteractor()->GetEventPosition()[0],
                this->GetInteractor()->GetEventPosition()[1],
                0, mView->GetRenderer() );
  picker->GetPickPosition( worldCoords );
  
  // Our world points are our RAS points.
  float RAS[3] = {0, 0, 0};
  for( int n = 0; n < 3; n++ )
    RAS[n] = worldCoords[n];
  
  // Did we hit a prop? Our only pickable prop is the surface.
  vtkProp3DCollection* props = picker->GetProp3Ds();
  if ( props->GetNumberOfItems() == 1 ) {
    
    props->InitTraversal();
    vtkProp* prop = props->GetNextProp();
    assert( prop == static_cast<vtkProp*>(mView->mSurfaceActor) );
    
    try {
      
      // Find the surface vertex at this point.
      float distance;
      int nVertex = mView->mCurrentSurfaceSource->
        FindVertexAtRAS( RAS, &distance );
      
      return nVertex;
    }
    catch(...) {
      // We clicked on the surface, but couldn't find a
      // vertex...
    }
  }

  // Didn't hit anything.
  return -1;
}

void
vtkKWQdecView::SetSurface ( vtkFSSurfaceSource* iSurface ) {

  //
  // Save the surface and make sure it's loaded, if it exists.
  //
  mCurrentSurfaceSource = iSurface;
  if( mCurrentSurfaceSource.GetPointer() )
    mCurrentSurfaceSource->Update();

  // Stick in our scalars now if we have them.
  if( mCurrentScalars.GetPointer() && mCurrentSurfaceSource.GetPointer() )
    mCurrentSurfaceSource->GetOutput()->GetPointData()->
      SetScalars( mCurrentScalars );

  //
  // Mapper for the surface.
  //
  mSurfaceMapper = NULL;
  if( mCurrentSurfaceSource.GetPointer() ) {
    mSurfaceMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mSurfaceMapper->SetInput( mCurrentSurfaceSource->GetOutput() );
  }

  // Set the color table if we have one.
  if( mCurrentScalarsColors.GetPointer() && mSurfaceMapper.GetPointer() )
    mSurfaceMapper->SetLookupTable( mCurrentScalarsColors );


  //
  // Surface actor.
  //
  mSurfaceActor = NULL;
  if( mSurfaceMapper.GetPointer() ) {
    mSurfaceActor = vtkSmartPointer<vtkActor>::New();
    mSurfaceActor->SetMapper( mSurfaceMapper );
    vtkSmartPointer<vtkProperty> property;
    property.TakeReference( mSurfaceActor->MakeProperty() );
    mSurfaceActor->SetProperty( property );
    property->SetColor( 0.6, 0.5, 0.5 );
    property.TakeReference( mSurfaceActor->MakeProperty() );
    mSurfaceActor->SetBackfaceProperty( property );
    property->BackfaceCullingOff();
  }

  //
  // The overlay surface is a surface with the same geometry but just
  // a little inflated so it sits over the main surface. This is the
  // inflator filter.
  //
  mSurfaceOverlayInflator = NULL;
  if( mCurrentSurfaceSource.GetPointer() ) {
    mSurfaceOverlayInflator = vtkSmartPointer<vtkInflatePolyData>::New();
    mSurfaceOverlayInflator->SetInflateFactor( 0.1 );
    mSurfaceOverlayInflator->
      SetInputConnection( mCurrentSurfaceSource->GetOutputPort() );
  }

  //
  // Mapper for the overlay surface.
  //
  mSurfaceOverlayMapper = NULL;
  if( mSurfaceOverlayInflator.GetPointer() ) {
    mSurfaceOverlayMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mSurfaceOverlayMapper->SetInput( mSurfaceOverlayInflator->GetOutput() );
  }
  
  //
  // Overlay surface actor.
  //
  mSurfaceOverlayActor = NULL;
  if( mSurfaceOverlayMapper.GetPointer() ) {
    mSurfaceOverlayActor = vtkSmartPointer<vtkActor>::New();
    mSurfaceOverlayActor->SetMapper( mSurfaceOverlayMapper );
    vtkSmartPointer<vtkProperty> property;
    property.TakeReference( mSurfaceOverlayActor->MakeProperty() );
    mSurfaceOverlayActor->SetProperty( property );
    property->SetOpacity( 0.0 );
    property->SetColor( 1, 1, 1 );
    mSurfaceOverlayActor->PickableOff();
  }

  // If we have surface overlay scalars and colors, set them now.
  if( mCurrentOverlayScalars.GetPointer() && 
      mCurrentOverlayColors.GetPointer() &&
      mSurfaceOverlayInflator.GetPointer() ) {

    mSurfaceOverlayInflator->Update(); // Update to generate the output
    mSurfaceOverlayInflator->GetOutput()->GetPointData()->
      SetScalars( mCurrentOverlayScalars );
    mSurfaceOverlayMapper->SetColorModeToMapScalars();
    mSurfaceOverlayMapper->UseLookupTableScalarRangeOn();
    mSurfaceOverlayMapper->SetLookupTable( mCurrentOverlayColors );
    mSurfaceOverlayActor->GetProperty()->SetOpacity( mOverlayOpacity );
  }    

  //
  // Remove current actors to clear the existing surface if
  // present. Then add the surface actors.
  //
  //this->GetRenderer()->SetBackground( 1, 1, 1 ); // white background
  this->GetRenderer()->RemoveAllViewProps();
  if( mSurfaceActor.GetPointer() )
    this->GetRenderer()->AddActor( mSurfaceActor );
  if( mSurfaceOverlayActor.GetPointer() )
    this->GetRenderer()->AddActor( mSurfaceOverlayActor );


  //
  // Cursor. Only create this if it doesn't already exist. Otherwise,
  // the existing one is fine, we'll just have to set the model
  // bounds.
  //
  mCursor = NULL;
  if( mCurrentSurfaceSource.GetPointer() ) {

    mCursor = vtkSmartPointer<vtkCursor3D>::New();
    mCursor->OutlineOff();
    mCursor->XShadowsOff();
    mCursor->YShadowsOff();
    mCursor->ZShadowsOff();
    
    // Cursor mapper.
    vtkSmartPointer<vtkPolyDataMapper> cursorMapper =
      vtkSmartPointer<vtkPolyDataMapper>::New();
    cursorMapper->SetInputConnection( mCursor->GetOutputPort() );
    
    // Cursor actor.
    vtkSmartPointer<vtkActor> cursorActor =
      vtkSmartPointer<vtkActor>::New();
    cursorActor->SetMapper( cursorMapper );
    vtkSmartPointer<vtkProperty> property;
    property.TakeReference( cursorActor->MakeProperty() );
    cursorActor->SetProperty( property );
    property->SetColor( 0, 1, 0 ); // green cross-hairs
    cursorActor->PickableOff();
    
    // Add the cursor actor.
    this->GetRenderer()->AddActor( cursorActor );
  
    // Set the bounds on the cursor.
    mCursor->SetModelBounds( mCurrentSurfaceSource->GetOutput()->GetBounds() );
  }

  //
  // Scalar bar actor.
  //
  if( mCurrentLegendColors.GetPointer() ) {

    // Create the scalar bar if not already create.
    this->InitializeScalarBar();

    // Populate the bar with the legend colors.
    if( mScalarBar.GetPointer() ) {
      mScalarBar->SetLookupTable( mCurrentLegendColors );
      
      // And the scalar bar actor.
      this->GetRenderer()->AddActor( mScalarBar );
    }
  }

  // If we had a cursor vertex, find the corresponding RAS point and
  // set the cursor axes.
  if( mCursor && -1 != mnCursorVertexIndex ) {
    float surfaceRAS[3];
    double surfaceRASd[3];
    mCurrentSurfaceSource->
      GetSurfaceRASAtVertex( mnCursorVertexIndex, surfaceRAS );
    for( int n = 0; n < 3; n++ )
      surfaceRASd[n] = surfaceRAS[n];
    mCursor->SetFocalPoint( surfaceRASd );
  }

  //
  // Selection points
  //
  // We use a list of points and lines as input to the poly data that
  // will draw the selection.
  // Make the polydata and set the points and lines.
  mSurfaceSelectionPolyData = vtkSmartPointer<vtkPolyData>::New();
  if( mSurfaceSelectionPolyData.GetPointer() ) {
    mSurfaceSelectionPolyData->SetPoints( mSurfaceSelectionPoints );
    mSurfaceSelectionPolyData->SetLines( mSurfaceSelectionLines );
  	
    // Mapper.
    vtkSmartPointer<vtkPolyDataMapper> pointsMapper = 
      vtkSmartPointer<vtkPolyDataMapper>::New();
    pointsMapper->SetInput( mSurfaceSelectionPolyData );
    
    // Actor.
    vtkSmartPointer<vtkActor> pointsActor =
      vtkSmartPointer<vtkActor>::New();
    pointsActor->SetMapper( pointsMapper );
    pointsActor->GetProperty()->SetColor( 0, 1, 0 );
    pointsActor->GetProperty()->SetPointSize( 2 );
    pointsActor->GetProperty()->SetLineWidth( 5 );
    pointsActor->PickableOff();
    
    // Add the selection actor.
    this->GetRenderer()->AddActor( pointsActor );
  }

  // Reset our path.
  this->ResetPath();

  // Render our new surface.
  this->Render();
}

void 
vtkKWQdecView::ResetView () {

  // Initial camera position to lateral view (assumes left hemi).
  this->GetRenderer()->GetActiveCamera()->SetPosition( -1, 0, 0 );
  this->GetRenderer()->GetActiveCamera()->SetFocalPoint( 0, 0, 0 );
  this->GetRenderer()->GetActiveCamera()->SetViewUp( 0, 0, 1 );

  // Reset the camera to show all the actors.
  this->GetRenderer()->ResetCamera();
  this->GetRenderer()->GetActiveCamera()->Zoom( mDefaultZoom );
  this->GetRenderer()->GetRenderWindow()->Render();

  // Grab the starting camera information so we can restore it later.
  this->GetRenderer()->GetActiveCamera()->GetPosition( mDefaultPosition );
  this->GetRenderer()->GetActiveCamera()->GetFocalPoint( mDefaultFocalPoint );
  this->GetRenderer()->GetActiveCamera()->GetViewUp( mDefaultViewUp );

}

void
vtkKWQdecView::RestoreView ( const char* isHemi ) {

  // Set the camera position to our original settings (assumes left hemi).
  this->GetRenderer()->GetActiveCamera()->SetPosition( mDefaultPosition );
  this->GetRenderer()->GetActiveCamera()->SetFocalPoint( mDefaultFocalPoint );
  this->GetRenderer()->GetActiveCamera()->SetViewUp( mDefaultViewUp );

  this->GetRenderer()->ResetCamera();
  this->GetRenderer()->GetActiveCamera()->Zoom( mDefaultZoom );
  this->GetRenderer()->GetRenderWindow()->Render();

  if( strcmp(isHemi, "rh") == 0 ) {
    this->AnimateCameraAzimuthNegative( );
    this->AnimateCameraAzimuthNegative( );
  }
}

void
vtkKWQdecView::ZoomBy ( float iFactor ) {
  this->GetRenderer()->GetActiveCamera()->Zoom( iFactor );
  this->GetRenderer()->GetRenderWindow()->Render();
}

bool
vtkKWQdecView::GetShowCursor () {
  if( mCursor ) 
    return mCursor->GetAxes();

  return false;
}
void
vtkKWQdecView::SetShowCursor ( bool ibShow ) {
  if( mCursor.GetPointer() ) {
    mCursor->SetAxes( ibShow );
    GetRenderer()->GetRenderWindow()->Render();
  }
}

void
vtkKWQdecView::SetSurfaceScalars ( vtkFloatArray* iScalars ) {

  // Save the pointer.
  mCurrentScalars = iScalars;

  // Set the scalars in the source.
  if( mCurrentSurfaceSource.GetPointer() )
    mCurrentSurfaceSource->GetOutput()->GetPointData()->
      SetScalars( mCurrentScalars );

}

void
vtkKWQdecView::SetSurfaceScalarsColors ( vtkScalarsToColors* iScalarsColors ) {

  // Save the pointer.
  mCurrentScalarsColors = iScalarsColors;

  // Set the color table in the surface mapper.
  if( mSurfaceMapper.GetPointer() )
    mSurfaceMapper->SetLookupTable( mCurrentScalarsColors );
}

void
vtkKWQdecView::SetSurfaceLookupScalars ( vtkFloatArray* iScalars ) {

  // Save the pointer.
  mCurrentLookupScalars = iScalars;
}

void
vtkKWQdecView::SetSurfaceOverlayScalarsAndColors 
( vtkFloatArray* iScalars,
  vtkScalarsToColors* iColors ) {

  // Save the pointers.
  mCurrentOverlayScalars = iScalars;
  mCurrentOverlayColors = iColors;

  // Set up the mapper and source.
  if( mSurfaceOverlayInflator.GetPointer() && 
      mCurrentOverlayScalars.GetPointer() ) {
    mSurfaceOverlayInflator->GetOutput()->GetPointData()->
      SetScalars( mCurrentOverlayScalars );
  }
  if( mSurfaceOverlayMapper.GetPointer() &&
      mCurrentOverlayColors.GetPointer() ) {
    mSurfaceOverlayMapper->SetLookupTable( mCurrentOverlayColors );
    mSurfaceOverlayMapper->UseLookupTableScalarRangeOn();
    mSurfaceOverlayActor->GetProperty()->SetOpacity( mOverlayOpacity );
  }
}

void
vtkKWQdecView::SetSurfaceOverlayOpacity ( double iOpacity ) {

  if( mSurfaceOverlayActor.GetPointer() ) {
    mOverlayOpacity = iOpacity;
    mSurfaceOverlayActor->GetProperty()->SetOpacity( mOverlayOpacity );
    this->Render();
  }
}

double
vtkKWQdecView::GetSurfaceOverlayOpacity () {

  return mOverlayOpacity;
}

void
vtkKWQdecView::SetSurfaceLegendColors ( vtkScalarsToColors* iLegendColors ) {

  // Save the pointer.
  mCurrentLegendColors = iLegendColors;

  // Set the color table in the bar. If we don't yet have a scalar
  // bar, make one and add it to the renderer.
  if( mCurrentLegendColors.GetPointer() ) {
    this->InitializeScalarBar();
    assert( mScalarBar );
    mScalarBar->SetLookupTable( mCurrentLegendColors );
  } else {
    if( mScalarBar.GetPointer() )
      this->DeleteScalarBar();
  }
}

void
vtkKWQdecView::SetROI ( vtkPolyData* iROIPolyData ) {

  // Save the pointer.
  mROIPolyData = iROIPolyData;

  // Initialize the actor.
  this->InitializeROI();

  // Render it.
  this->Render();
}

void
vtkKWQdecView::SetAnnotationMessage ( const char* isMessage ) {

  // Set the message in the top left corner.
  this->GetCornerAnnotation()->SetText( 2, isMessage );
  this->Render();
}

void
vtkKWQdecView::SetShowLegend ( int ibShow ) {

  if( mScalarBar.GetPointer() ) {
    mScalarBar->SetVisibility( ibShow );
  }
}

void
vtkKWQdecView::SetShowAnnotation ( int ibShow ) {

  this->SetAnnotationsVisibility( ibShow );
}

void
vtkKWQdecView::SetVertexAnnotationLookup 
( QdecVertexAnnotationLookup* iLookup ) {

  mAnnotationLookup = iLookup;
}

vtkPoints*
vtkKWQdecView::GetSurfaceSelectionPoints () {

  assert( mSurfaceSelectionPoints.GetPointer() );

  return mSurfaceSelectionPoints;
}

void
vtkKWQdecView::SelectSurfaceVertex ( int inVertex ) {

  if( NULL == mCurrentSurfaceSource.GetPointer() )
    throw runtime_error( "No surface loaded." );

  if( inVertex >= mCurrentSurfaceSource->GetNumberOfVertices() )
    throw runtime_error( "Invalid vertex number." );

  // If this is -1, just unselect whatever is current.
  if( -1 == inVertex ) {

    // Clear the annotation.
    this->GetCornerAnnotation()->SetText( 0, "" );

    // Clear the GDF if loaded.
    if( mVertexPlot && mVertexPlot->IsLoaded() ) {
      mVertexPlot->BeginPointList();
      mVertexPlot->EndPointList();
    }

  } else {

    // Start building a label for the annotation.
    stringstream ssLabel;
    
    // Put the location in the label.
    float surfaceRAS[3];
    mCurrentSurfaceSource->GetSurfaceRASAtVertex( inVertex, surfaceRAS );
    ssLabel << "(" << fixed << setprecision(2)
            << surfaceRAS[0] << ", "
            << surfaceRAS[1] << ", " 
            << surfaceRAS[2] << ")";
    
    // Add the info to the label.
    ssLabel << " Vertex #" << inVertex;
    
    // If we have lookup scalars, print the value.
    if( mCurrentLookupScalars ) {
      ssLabel << " value " << mCurrentLookupScalars->GetTuple1( inVertex );
    }
    
    // If we have an annotatoin lookup, get the string and add it to
    // the label.
    if( mAnnotationLookup ) {
      ssLabel << " " << mAnnotationLookup->GetAnnotationForVertex(inVertex);
    }
    
    // If our GDF is loaded, plot this point and set the info.
    if( mVertexPlot && mVertexPlot->IsLoaded() ) {
      mVertexPlot->SetPoint( inVertex );
      mVertexPlot->SetInfo( ssLabel.str().c_str() );
    }
  
    // Set the label in the lower left annotation.
    this->GetCornerAnnotation()->SetText( 0, ssLabel.str().c_str() );
    
    // Set the cursor focal point.
    if( mCursor.GetPointer() ) {
      double worldCoords[3];
      for( int n = 0; n < 3; n ++ )
        worldCoords[n] = static_cast<double>( surfaceRAS[n] );
      mCursor->SetFocalPoint( worldCoords );
    }
  }

  // Save the index of this vertex.
  mnCursorVertexIndex = inVertex;

  // Broadcast our event.
  SurfaceVertexInformation info;
  info.mnVertexIndex = mnCursorVertexIndex;
  this->InvokeEvent( QdecEvents::UserSelectedVertex, 
                     static_cast<void*>(&info) );
    
  // Render the view.
  this->Render();
}

void
vtkKWQdecView::InitializeScalarBar () {

  // If we already have one, delete the actor from the view and delete
  // the bar.
  if( mScalarBar.GetPointer() )
    this->DeleteScalarBar();

  // Make a new one.
  mScalarBar = vtkSmartPointer<vtkScalarBarActor>::New();
#if 0
  // vertical, on right
  mScalarBar->SetOrientationToVertical();
  mScalarBar->GetLabelTextProperty()->SetFontSize( 5 );
  mScalarBar->SetPosition( 0.9, 0.1 );
  mScalarBar->SetWidth( 0.1 );
#else
  // horizonal, at bottom right
  mScalarBar->SetOrientationToHorizontal();
  mScalarBar->GetLabelTextProperty()->SetFontSize( 4 );
  mScalarBar->SetPosition( 0.66, 0.0 );
  mScalarBar->SetHeight( 0.1 );
  mScalarBar->SetWidth( 0.3 );
#endif
  mScalarBar->PickableOff();

  // Add it to the view.
  this->GetRenderer()->AddActor( mScalarBar );
}

void
vtkKWQdecView::DeleteScalarBar () {

  assert( mScalarBar.GetPointer() );

  if( this->GetRenderer()->HasViewProp( mScalarBar ) )
    this->GetRenderer()->RemoveViewProp( mScalarBar );
  mScalarBar = NULL;
}

void
vtkKWQdecView::InitializeROI () {

  // Create our inflator. This makes sure the ROI is drawn just above
  // the surface.
  vtkSmartPointer<vtkInflatePolyData> inflator =
    vtkSmartPointer<vtkInflatePolyData>::New();
  inflator->SetInput( mROIPolyData );
  inflator->SetInflateFactor( 0.1 );

  // Make our mapper.
  mROIMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  mROIMapper->SetInput( inflator->GetOutput() );
  
  // Make our actor.
  if( mROIActor.GetPointer() )
    this->GetRenderer()->RemoveActor( mROIActor );
  mROIActor = vtkActor::New();
  mROIActor->SetMapper( mROIMapper );
  mROIActor->PickableOff();
  vtkSmartPointer<vtkProperty> property;
  property.TakeReference( mROIActor->MakeProperty() );
  mROIActor->SetProperty( property );
  property->SetColor( 0.5, 0.5, 0.9 );

  // Add the actor to the view.
  this->GetRenderer()->AddActor( mROIActor );
}    

void 
vtkKWQdecView::AnimateCameraElevatePositive () {

  double rate = this->GetRenderWindowInteractor()->GetDesiredUpdateRate();
  this->GetRenderWindowInteractor()->
    SetDesiredUpdateRate( kcAnimateSteps / kAnimateTimeInSeconds );
  for( double nStep = 1.0; nStep <= kcAnimateSteps; nStep += 1.0 ) {
    this->GetRenderer()->GetActiveCamera()->
      Elevation( kAnimationDegrees/kcAnimateSteps );
    this->GetRenderer()->GetActiveCamera()->OrthogonalizeViewUp();
    this->Render();
  }
  this->GetRenderWindowInteractor()->SetDesiredUpdateRate( rate );
}

void 
vtkKWQdecView::AnimateCameraElevateNegative () {

  double rate = this->GetRenderWindowInteractor()->GetDesiredUpdateRate();
  this->GetRenderWindowInteractor()->
    SetDesiredUpdateRate( kcAnimateSteps / kAnimateTimeInSeconds );
  for( double nStep = 1.0; nStep <= kcAnimateSteps; nStep += 1.0 ) {
    this->GetRenderer()->GetActiveCamera()->
      Elevation( -kAnimationDegrees/kcAnimateSteps );
    this->GetRenderer()->GetActiveCamera()->OrthogonalizeViewUp();
    this->Render();
  }
  this->GetRenderWindowInteractor()->SetDesiredUpdateRate( rate );
}

void 
vtkKWQdecView::AnimateCameraAzimuthNegative () {

  double rate = this->GetRenderWindowInteractor()->GetDesiredUpdateRate();
  this->GetRenderWindowInteractor()->
    SetDesiredUpdateRate( kcAnimateSteps / kAnimateTimeInSeconds );
  for( double nStep = 1.0; nStep <= kcAnimateSteps; nStep += 1.0 ) {
    this->GetRenderer()->GetActiveCamera()->
      Azimuth( -kAnimationDegrees/kcAnimateSteps );
    this->GetRenderer()->GetActiveCamera()->OrthogonalizeViewUp();
    this->Render();
  }
  this->GetRenderWindowInteractor()->SetDesiredUpdateRate( rate );
}

void 
vtkKWQdecView::AnimateCameraAzimuthPositive () {

  double rate = this->GetRenderWindowInteractor()->GetDesiredUpdateRate();
  this->GetRenderWindowInteractor()->
    SetDesiredUpdateRate( kcAnimateSteps / kAnimateTimeInSeconds );
  for( double nStep = 1.0; nStep <= kcAnimateSteps; nStep += 1.0 ) {
    this->GetRenderer()->GetActiveCamera()->
      Azimuth( kAnimationDegrees/kcAnimateSteps );
    this->GetRenderer()->GetActiveCamera()->OrthogonalizeViewUp();
    this->Render();
  }
  this->GetRenderWindowInteractor()->SetDesiredUpdateRate( rate );
}

void 
vtkKWQdecView::AnimateCameraRollNegative () {

  double rate = this->GetRenderWindowInteractor()->GetDesiredUpdateRate();
  this->GetRenderWindowInteractor()->
    SetDesiredUpdateRate( kcAnimateSteps / kAnimateTimeInSeconds );
  for( double nStep = 1.0; nStep <= kcAnimateSteps; nStep += 1.0 ) {
    this->GetRenderer()->GetActiveCamera()->
      Roll( -kAnimationDegrees/kcAnimateSteps );
    this->Render();
  }
  this->GetRenderWindowInteractor()->SetDesiredUpdateRate( rate );
}

void 
vtkKWQdecView::AnimateCameraRollPositive () {

  double rate = this->GetRenderWindowInteractor()->GetDesiredUpdateRate();
  this->GetRenderWindowInteractor()->
    SetDesiredUpdateRate( kcAnimateSteps / kAnimateTimeInSeconds );
  for( double nStep = 1.0; nStep <= kcAnimateSteps; nStep += 1.0 ) {
    this->GetRenderer()->GetActiveCamera()->
      Roll( kAnimationDegrees/kcAnimateSteps );
    this->Render();
  }
  this->GetRenderWindowInteractor()->SetDesiredUpdateRate( rate );
}

void
vtkKWQdecView::ResetPath () {

  assert( mSurfaceSelectionPoints.GetPointer() );
  assert( mSurfaceSelectionLines.GetPointer() );

  // Clear the lines.
  mSurfaceSelectionPoints->Reset();
  mSurfaceSelectionLines->Reset();
}

void
vtkKWQdecView::AddVertexToPath ( int inVertex ) {

  assert( mSurfaceSelectionPoints.GetPointer() );
  assert( mCurrentSurfaceSource.GetPointer() );

  // Get a point from the vertex index.
  double surfaceRAS[3];
  mCurrentSurfaceSource->GetSurfaceRASAtVertex( inVertex, surfaceRAS );

  // Just insert this point.
  mSurfaceSelectionPoints->InsertNextPoint( surfaceRAS );
}

void
vtkKWQdecView::RebuildPathLine () {

  assert( mSurfaceSelectionLines.GetPointer() );
  assert( mSurfaceSelectionPoints.GetPointer() );

  // Reset the lines first.
  mSurfaceSelectionLines->Reset();

  // For each point in our list of points, add a line segment for
  // those that point and the next one.
  for( int n = 0; n+1 < mSurfaceSelectionPoints->GetNumberOfPoints(); n++ ) {
    vtkIdType pts[2];
    pts[0] = n;
    pts[1] = n+1;
    mSurfaceSelectionLines->InsertNextCell( 2, pts );
  }

  // Let the poly data know we've changed.
  mSurfaceSelectionPolyData->Modified();
  mSurfaceSelectionPolyData->Update();
}
