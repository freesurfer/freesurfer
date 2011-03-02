/**
 * @file  vtkOrientMRIInteractorStyleView2D.cxx
 * @brief Interactor style for the 2D view
 *
 * This interactor style allows the camera to pan and dolly so the
 * user can see the plane better, but doesn't allow him to rotate the
 * plane. Instead, the left button does a window/level operation. So
 * this is just a subclass of vtkInteractorStyleTrackballCamera with
 * the left button overridden. The code for window/level interaction
 * comes from vtkImagePlaneWidget.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:35 $
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

#include <assert.h>

#include "vtkOrientMRIInteractorStyleView2D.h"

#include "OrientMRIEvents.h"
#include "vtkAbstractPicker.h"
#include "vtkCommand.h"
#include "vtkCallbackCommand.h"
#include "vtkMatrix4x4.h"
#include "vtkObjectFactory.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkWindowLevelLookupTable.h"

vtkCxxRevisionMacro(vtkOrientMRIInteractorStyleView2D, "$Revision: 1.5 $");
vtkStandardNewMacro(vtkOrientMRIInteractorStyleView2D);

vtkOrientMRIInteractorStyleView2D::vtkOrientMRIInteractorStyleView2D () :
  OriginalWindow( 1.0 ),
  ActiveOrthoLine( -1 ) {
}

vtkOrientMRIInteractorStyleView2D::~vtkOrientMRIInteractorStyleView2D () {
}

void
vtkOrientMRIInteractorStyleView2D::SetWindowLevelTable ( vtkWindowLevelLookupTable* iTable ) {

  // Save a pointer.
  this->LookupTable = iTable;

  // If we got a table, get the range. It will be used as a scaling
  // factor.
  if( this->LookupTable.GetPointer() ) {
    this->OriginalWindow = this->LookupTable->GetWindow();
  }
}

void
vtkOrientMRIInteractorStyleView2D::
SetPointsAndLines ( double* iPointsButton1Start, 
		    double* iPointsButton1End,
		    double* iPointsButton2Start, 
		    double* iPointsButton2End,
		    double* iPointsButton3Start, 
		    double* iPointsButton3End ) {

  this->Points[0][0] = iPointsButton1Start;
  this->Points[0][1] = iPointsButton1End;
  this->Points[1][0] = iPointsButton2Start;
  this->Points[1][1] = iPointsButton2End;
  this->Points[2][0] = iPointsButton3Start;
  this->Points[2][1] = iPointsButton3End;
}

void
vtkOrientMRIInteractorStyleView2D::OnMouseMove () {

  // If we're in the middle of a window/level or ortho line edit
  // operation, call our function to do that, otherwise just let the
  // superclass handle it.
  if( IS_WindowLevel == this->State ) {

    this->WindowLevel();
    this->InvokeEvent( vtkCommand::InteractionEvent, NULL );
    
  } else if( IS_OrthoLineEdit == this->State ) {

    this->EditOrthoLine();
    this->InvokeEvent( vtkCommand::InteractionEvent, NULL );
    
  } else {

    return this->vtkInteractorStyleTrackballCamera::OnMouseMove();
  }
}

void
vtkOrientMRIInteractorStyleView2D::OnLeftButtonDown () {
  int x = this->Interactor->GetEventPosition()[0];
  int y = this->Interactor->GetEventPosition()[1];

  // Find a renderer.
  this->FindPokedRenderer( x, y );
  if( NULL == this->CurrentRenderer )
    return;

  // If this is a control click, start doing a line ortho
  // edit. Otherwise, do a window/level operation.
  if( this->Interactor->GetControlKey() ) {
#if ((VTK_MAJOR_VERSION == 5)&&(VTK_MINOR_VERSION > 0))
    this->GrabFocus( this->EventCallbackCommand );
#endif
    this->StartOrthoLineEdit( 0 );
  } else {
#if ((VTK_MAJOR_VERSION == 5)&&(VTK_MINOR_VERSION > 0))
    this->GrabFocus( this->EventCallbackCommand );
#endif
    this->StartWindowLevel();
  }
}

void
vtkOrientMRIInteractorStyleView2D::OnLeftButtonUp () {
  
  // If we're in the middle of a window/level or ortho line
  // edit operation, end it. Otherwise just let the superclass handle it.
  if( IS_WindowLevel == this->State ) {
    
    this->EndWindowLevel();
#if ((VTK_MAJOR_VERSION == 5)&&(VTK_MINOR_VERSION > 0))
    if( this->Interactor )
      this->ReleaseFocus();
#endif

  } else if( IS_OrthoLineEdit == this->State ) {
    
    this->EndOrthoLineEdit();
#if ((VTK_MAJOR_VERSION == 5)&&(VTK_MINOR_VERSION > 0))
    if( this->Interactor )
      this->ReleaseFocus();
#endif

  } else {
    
    return this->vtkInteractorStyleTrackballCamera::OnLeftButtonUp();
  }
}

void
vtkOrientMRIInteractorStyleView2D::OnMiddleButtonDown () {
  int x = this->Interactor->GetEventPosition()[0];
  int y = this->Interactor->GetEventPosition()[1];

  // Find a renderer.
  this->FindPokedRenderer( x, y );
  if( NULL == this->CurrentRenderer )
    return;

  // If this is a control click, start doing a line ortho
  // edit. 
  if( this->Interactor->GetControlKey() ) {
#if ((VTK_MAJOR_VERSION == 5)&&(VTK_MINOR_VERSION > 0))
    this->GrabFocus( this->EventCallbackCommand );
#endif
    this->StartOrthoLineEdit( 1 );
  } else {
    return this->vtkInteractorStyleTrackballCamera::OnMiddleButtonDown();
  }
}

void
vtkOrientMRIInteractorStyleView2D::OnMiddleButtonUp () {
  
  // If we're in the middle of a ortho line edit operation, end
  // it. Otherwise just let the superclass handle it.
  if( IS_OrthoLineEdit == this->State ) {
    
    this->EndOrthoLineEdit();
#if ((VTK_MAJOR_VERSION == 5)&&(VTK_MINOR_VERSION > 0))
    if( this->Interactor )
      this->ReleaseFocus();
#endif

  } else {
    
    return this->vtkInteractorStyleTrackballCamera::OnMiddleButtonUp();
  }
}


void
vtkOrientMRIInteractorStyleView2D::OnRightButtonDown () {
  int x = this->Interactor->GetEventPosition()[0];
  int y = this->Interactor->GetEventPosition()[1];

  // Find a renderer.
  this->FindPokedRenderer( x, y );
  if( NULL == this->CurrentRenderer )
    return;

  // If this is a control click, start doing a line ortho
  // edit. 
  if( this->Interactor->GetControlKey() ) {
#if ((VTK_MAJOR_VERSION == 5)&&(VTK_MINOR_VERSION > 0))
    this->GrabFocus( this->EventCallbackCommand );
#endif
    this->StartOrthoLineEdit( 2 );
  }  else {
    return this->vtkInteractorStyleTrackballCamera::OnRightButtonDown();
  }
}

void
vtkOrientMRIInteractorStyleView2D::OnRightButtonUp () {
  
  // If we're in the middle of a ortho line edit operation, end
  // it. Otherwise just let the superclass handle it.
  if( IS_OrthoLineEdit == this->State ) {
    
    this->EndOrthoLineEdit();
#if ((VTK_MAJOR_VERSION == 5)&&(VTK_MINOR_VERSION > 0))
    if( this->Interactor )
      this->ReleaseFocus();
#endif

  } else {
    
    return this->vtkInteractorStyleTrackballCamera::OnRightButtonUp();
  }
}

void
vtkOrientMRIInteractorStyleView2D::WindowLevel () {
  
  if( NULL == this->LookupTable.GetPointer() )
    return;

  int x = this->Interactor->GetEventPosition()[0];
  int y = this->Interactor->GetEventPosition()[1];
  
  // Get the current window and level.
  double window = this->LookupTable->GetWindow();
  double level = this->LookupTable->GetLevel();
  
  // The original window acts as a factor for adjusting the current
  // values.
  double owin = this->OriginalWindow;

  // Adjust the window and level using mouse movements (level on x,
  // window on y), and apply the original window factor.
  level =
    level + (x - this->Interactor->GetLastEventPosition()[0])* owin/500.0;
  window =
     window + (this->Interactor->GetLastEventPosition()[1] - y)* owin/250.0;
  
  // Set the new values and redraw.
  this->LookupTable->SetWindow( window );
  this->LookupTable->SetLevel( level );
  
  this->Interactor->Render();
}

void
vtkOrientMRIInteractorStyleView2D::EditOrthoLine () {

  assert( this->Points[this->ActiveOrthoLine][0] );
  assert( this->Points[this->ActiveOrthoLine][1] );
  
  // This is the current position.
  int x = this->Interactor->GetEventPosition()[0];
  int y = this->Interactor->GetEventPosition()[1];

  // Pick a world point.
  vtkAbstractPicker* picker = this->Interactor->GetPicker();
  assert( picker );
  double currentWorldCoords[3] = {0,0,0};
  picker->Pick( x, y, 0, this->GetCurrentRenderer() );
  picker->GetPickPosition( currentWorldCoords );

  // Copy our points into the pointers we got.
  this->Points[this->ActiveOrthoLine][0][0] = this->ActiveOrthoLineStart[0];
  this->Points[this->ActiveOrthoLine][0][1] = this->ActiveOrthoLineStart[1];
  this->Points[this->ActiveOrthoLine][0][2] = this->ActiveOrthoLineStart[2];
  this->Points[this->ActiveOrthoLine][1][0] = currentWorldCoords[0];
  this->Points[this->ActiveOrthoLine][1][1] = currentWorldCoords[1];
  this->Points[this->ActiveOrthoLine][1][2] = currentWorldCoords[2];
 
  // Broadcast our event.
  this->InvokeEvent( OrientMRIEvents::OrthoLineChanged, 
		     static_cast<void*>(&this->ActiveOrthoLine) );

  this->Interactor->Render();

}

void
vtkOrientMRIInteractorStyleView2D::StartWindowLevel () {

  if( this->State != VTKIS_NONE )
    return;

  this->StartState( IS_WindowLevel );
}

void
vtkOrientMRIInteractorStyleView2D::EndWindowLevel () {

  if( this->State != IS_WindowLevel )
    return;

  this->StopState();
}

void
vtkOrientMRIInteractorStyleView2D::StartOrthoLineEdit ( int iWhich ) {

  if( this->State != VTKIS_NONE )
    return;

  // Remember which line we're editing.
  this->ActiveOrthoLine = iWhich;

  // Remember the mouse down position.
  this->ActiveOrthoLineStart[0] = this->Interactor->GetEventPosition()[0];
  this->ActiveOrthoLineStart[1] = this->Interactor->GetEventPosition()[1];

  // Get the world coords here.
  vtkAbstractPicker* picker = this->Interactor->GetPicker();
  assert( picker );
  picker->Pick( this->ActiveOrthoLineStart[0],
		this->ActiveOrthoLineStart[1], 0, 
		this->GetCurrentRenderer() );
  picker->GetPickPosition( this->ActiveOrthoLineStart );

  this->StartState( IS_OrthoLineEdit );
}

void
vtkOrientMRIInteractorStyleView2D::EndOrthoLineEdit () {

  if( this->State != IS_OrthoLineEdit )
    return;

  this->StopState();
}
