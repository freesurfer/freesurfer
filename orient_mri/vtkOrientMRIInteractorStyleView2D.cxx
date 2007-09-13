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
 *    $Author: kteich $
 *    $Date: 2007/09/13 20:58:21 $
 *    $Revision: 1.1 $
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

#include "vtkOrientMRIInteractorStyleView2D.h"

#include "vtkCommand.h"
#include "vtkCallbackCommand.h"
#include "vtkMatrix4x4.h"
#include "vtkObjectFactory.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkWindowLevelLookupTable.h"

vtkCxxRevisionMacro(vtkOrientMRIInteractorStyleView2D, "$Revision: 1.1 $");
vtkStandardNewMacro(vtkOrientMRIInteractorStyleView2D);

vtkOrientMRIInteractorStyleView2D::vtkOrientMRIInteractorStyleView2D () :
  OriginalWindow( 1.0 ) {
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
vtkOrientMRIInteractorStyleView2D::OnMouseMove () {

  int x = this->Interactor->GetEventPosition()[0];
  int y = this->Interactor->GetEventPosition()[1];

  // If we're in the middle of a window/level operation, call our
  // function to do that, otherwise just let the superclass handle it.
  if( IS_WindowLevel == this->State ) {

    this->FindPokedRenderer( x, y );
    this->WindowLevel();
    this->InvokeEvent( vtkCommand::InteractionEvent, NULL );
    
  } else {
    
    return this->vtkInteractorStyleTrackballCamera::OnMouseMove();
  }
}

void
vtkOrientMRIInteractorStyleView2D::OnLeftButtonDown () {
  int x = this->Interactor->GetEventPosition()[0];
  int y = this->Interactor->GetEventPosition()[1];

  // If we have a renderer, grab focus and start the window/level
  // operation.
  this->FindPokedRenderer( x, y );
  if( NULL == this->CurrentRenderer )
    return;

  this->GrabFocus( this->EventCallbackCommand );
  this->StartWindowLevel();
}

void
vtkOrientMRIInteractorStyleView2D::OnLeftButtonUp () {
  
  // If we're in the middle of a window/level operation, end
  // it. Otherwise just let the superclass handle it.
  if( IS_WindowLevel == this->State ) {
    
    this->EndWindowLevel();
    if( this->Interactor )
      this->ReleaseFocus();

  } else {
    
    return this->vtkInteractorStyleTrackballCamera::OnMouseMove();
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
