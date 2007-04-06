/**
 * @file  vtkScubaInteractorStyle.cxx
 * @brief vtk interactor style to hook into our tool class.
 *
 *  This is an implementation of the VTK interactor style method, a
 * listener that responds to interaction events. This class takes
 * those events and passes them to the vtkKWScubaView which then
 * passes them to a tool. Also tells the window to select the view
 * that was clicked.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: kteich $
 *    $Date: 2007/04/06 22:23:06 $
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


#include "vtkScubaInteractorStyle.h"
#include "vtkObjectFactory.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkRenderWindow.h"
#include "vtkAbstractPicker.h"
#include "vtkRenderer.h"
#include "vtkKWScubaWindow.h"
#include "vtkKWScubaView.h"

vtkStandardNewMacro( vtkScubaInteractorStyle );
vtkCxxRevisionMacro( vtkScubaInteractorStyle, "$Revision: 1.1 $" );

vtkScubaInteractorStyle::vtkScubaInteractorStyle () :
    mWindow(NULL) {}

vtkScubaInteractorStyle::~vtkScubaInteractorStyle () {}

void
vtkScubaInteractorStyle::SetWindow ( vtkKWScubaWindow* iWindow ) {

  mWindow = iWindow;
}

void
vtkScubaInteractorStyle::OnMouseMove () {

  this->FindPokedRenderer( this->Interactor->GetEventPosition()[0],
                           this->Interactor->GetEventPosition()[1] );
  if ( NULL == this->CurrentRenderer ) {
    return;
  }

  // See if we can find a view for this event.
  vtkRenderWindow* currentWindow = this->CurrentRenderer->GetRenderWindow();
  vtkKWScubaView* view =
    vtkKWScubaView::GetViewFromRenderWindow( currentWindow );
  if ( NULL == view ) {
    return;
  }

  // Notify the view of the event.
  view->MouseMoveEvent( mWindow, this->Interactor->GetEventPosition()  );

  // After the event is handled, notify the main window.
  mWindow->EventDone();
}

void
vtkScubaInteractorStyle::OnLeftButtonDown () {

  this->FindPokedRenderer( this->Interactor->GetEventPosition()[0],
                           this->Interactor->GetEventPosition()[1] );
  if ( NULL == this->CurrentRenderer ) {
    return;
  }

  // See if we can find a view for this event.
  vtkRenderWindow* currentWindow = this->CurrentRenderer->GetRenderWindow();
  vtkKWScubaView* view =
    vtkKWScubaView::GetViewFromRenderWindow( currentWindow );
  if ( NULL == view ) {
    return;
  }
  
  // Set the view clicked.
  mWindow->SetCurrentView( *view );

  // Notify the view of the event.
  view->LeftButtonDownEvent( mWindow, this->Interactor->GetEventPosition()  );

  // After the event is handled, notify the main window.
  mWindow->EventDone();
}

void
vtkScubaInteractorStyle::OnLeftButtonUp () {

  this->FindPokedRenderer( this->Interactor->GetEventPosition()[0],
                           this->Interactor->GetEventPosition()[1] );
  if ( NULL == this->CurrentRenderer ) {
    return;
  }

  // See if we can find a view for this event.
  vtkRenderWindow* currentWindow = this->CurrentRenderer->GetRenderWindow();
  vtkKWScubaView* view =
    vtkKWScubaView::GetViewFromRenderWindow( currentWindow );
  if ( NULL == view ) {
    return;
  }

  // Notify the view of the event.
  view->LeftButtonUpEvent( mWindow, this->Interactor->GetEventPosition()  );

  // After the event is handled, notify the main window.
  mWindow->EventDone();
}

void
vtkScubaInteractorStyle::OnMiddleButtonDown () {

  this->FindPokedRenderer( this->Interactor->GetEventPosition()[0],
                           this->Interactor->GetEventPosition()[1] );
  if ( NULL == this->CurrentRenderer ) {
    return;
  }

  // See if we can find a view for this event.
  vtkRenderWindow* currentWindow = this->CurrentRenderer->GetRenderWindow();
  vtkKWScubaView* view =
    vtkKWScubaView::GetViewFromRenderWindow( currentWindow );
  if ( NULL == view ) {
    return;
  }

  // Set the view clicked.
  mWindow->SetCurrentView( *view );

  // Notify the view of the event.
  view->MiddleButtonDownEvent( mWindow, this->Interactor->GetEventPosition()  );

  // After the event is handled, notify the main window.
  mWindow->EventDone();
}

void
vtkScubaInteractorStyle::OnMiddleButtonUp () {

  this->FindPokedRenderer( this->Interactor->GetEventPosition()[0],
                           this->Interactor->GetEventPosition()[1] );
  if ( NULL == this->CurrentRenderer ) {
    return;
  }

  // See if we can find a view for this event.
  vtkRenderWindow* currentWindow = this->CurrentRenderer->GetRenderWindow();
  vtkKWScubaView* view =
    vtkKWScubaView::GetViewFromRenderWindow( currentWindow );
  if ( NULL == view ) {
    return;
  }

  // Notify the view of the event.
  view->MiddleButtonUpEvent( mWindow, this->Interactor->GetEventPosition()  );

  // After the event is handled, notify the main window.
  mWindow->EventDone();
}

void
vtkScubaInteractorStyle::OnRightButtonDown () {

  this->FindPokedRenderer( this->Interactor->GetEventPosition()[0],
                           this->Interactor->GetEventPosition()[1] );
  if ( NULL == this->CurrentRenderer ) {
    return;
  }

  // See if we can find a view for this event.
  vtkRenderWindow* currentWindow = this->CurrentRenderer->GetRenderWindow();
  vtkKWScubaView* view =
    vtkKWScubaView::GetViewFromRenderWindow( currentWindow );
  if ( NULL == view ) {
    return;
  }

  // Set the view clicked.
  mWindow->SetCurrentView( *view );

  // Notify the view of the event.
  view->RightButtonDownEvent( mWindow, this->Interactor->GetEventPosition()  );

  // After the event is handled, notify the main window.
  mWindow->EventDone();
}

void
vtkScubaInteractorStyle::OnRightButtonUp () {

  this->FindPokedRenderer( this->Interactor->GetEventPosition()[0],
                           this->Interactor->GetEventPosition()[1] );
  if ( NULL == this->CurrentRenderer ) {
    return;
  }

  // See if we can find a view for this event.
  vtkRenderWindow* currentWindow = this->CurrentRenderer->GetRenderWindow();
  vtkKWScubaView* view =
    vtkKWScubaView::GetViewFromRenderWindow( currentWindow );
  if ( NULL == view ) {
    return;
  }

  // Notify the view of the event.
  view->RightButtonUpEvent( mWindow, this->Interactor->GetEventPosition()  );

  // After the event is handled, notify the main window.
  mWindow->EventDone();
}

void
vtkScubaInteractorStyle::OnKeyDown () {

  this->FindPokedRenderer( this->Interactor->GetEventPosition()[0],
                           this->Interactor->GetEventPosition()[1] );
  if ( NULL == this->CurrentRenderer ) {
    return;
  }

  // See if we can find a view for this event.
  vtkRenderWindow* currentWindow = this->CurrentRenderer->GetRenderWindow();
  vtkKWScubaView* view =
    vtkKWScubaView::GetViewFromRenderWindow( currentWindow );
  if ( NULL == view ) {
    return;
  }

  // Set the view clicked.
  mWindow->SetCurrentView( *view );

  // Notify the view of the event.
  view->KeyDownEvent( mWindow, this->Interactor->GetEventPosition()  );

  // After the event is handled, notify the main window.
  mWindow->EventDone();
}

void
vtkScubaInteractorStyle::OnKeyUp () {

  this->FindPokedRenderer( this->Interactor->GetEventPosition()[0],
                           this->Interactor->GetEventPosition()[1] );
  if ( NULL == this->CurrentRenderer ) {
    return;
  }

  // See if we can find a view for this event.
  vtkRenderWindow* currentWindow = this->CurrentRenderer->GetRenderWindow();
  vtkKWScubaView* view =
    vtkKWScubaView::GetViewFromRenderWindow( currentWindow );
  if ( NULL == view ) {
    return;
  }

  // Notify the view of the event.
  view->KeyUpEvent( mWindow, this->Interactor->GetEventPosition()  );

  // After the event is handled, notify the main window.
  mWindow->EventDone();
}

void
vtkScubaInteractorStyle::OnEnter () {

  this->FindPokedRenderer( this->Interactor->GetEventPosition()[0],
                           this->Interactor->GetEventPosition()[1] );
  if ( NULL == this->CurrentRenderer ) {
    return;
  }

  // See if we can find a view for this event.
  vtkRenderWindow* currentWindow = this->CurrentRenderer->GetRenderWindow();
  vtkKWScubaView* view =
    vtkKWScubaView::GetViewFromRenderWindow( currentWindow );
  if ( NULL == view ) {
    return;
  }

  // Notify the view of the event.
  view->EnterEvent( mWindow, this->Interactor->GetEventPosition()  );

  // After the event is handled, notify the main window.
  mWindow->EventDone();
}

void
vtkScubaInteractorStyle::OnLeave () {

  this->FindPokedRenderer( this->Interactor->GetEventPosition()[0],
                           this->Interactor->GetEventPosition()[1] );
  if ( NULL == this->CurrentRenderer ) {
    return;
  }

  // See if we can find a view for this event.
  vtkRenderWindow* currentWindow = this->CurrentRenderer->GetRenderWindow();
  vtkKWScubaView* view =
    vtkKWScubaView::GetViewFromRenderWindow( currentWindow );
  if ( NULL == view ) {
    return;
  }

  // Notify the view of the event.
  view->LeaveEvent( mWindow, this->Interactor->GetEventPosition()  );

  // After the event is handled, notify the main window.
  mWindow->EventDone();
}

