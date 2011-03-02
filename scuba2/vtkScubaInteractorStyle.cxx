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
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:40 $
 *    $Revision: 1.2 $
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


#include "vtkScubaInteractorStyle.h"
#include "vtkObjectFactory.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkRenderWindow.h"
#include "vtkAbstractPicker.h"
#include "vtkRenderer.h"
#include "vtkKWScubaWindow.h"
#include "vtkKWScubaView.h"

vtkStandardNewMacro( vtkScubaInteractorStyle );
vtkCxxRevisionMacro( vtkScubaInteractorStyle, "$Revision: 1.2 $" );

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

