/**
 * @file  vtkKWScubaToolNavigate.cxx
 * @brief Navigate tool
 *
 * Implements panning, rotating, 2DRASZ shifts, and zooming. Works in
 * 2D and 3D.
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


#include "vtkKWScubaToolNavigate.h"
#include "vtkObjectFactory.h"
#include "vtkKWScubaView.h"
#include "vtkRenderWindow.h"
#include "vtkCommand.h"
#include "vtkRenderWindowInteractor.h"

using namespace std;

vtkStandardNewMacro( vtkKWScubaToolNavigate );
vtkCxxRevisionMacro( vtkKWScubaToolNavigate, "$Revision: 1.2 $" );

vtkKWScubaToolNavigate::vtkKWScubaToolNavigate () {

  msLabel = "Navigate";
}

vtkKWScubaToolNavigate::~vtkKWScubaToolNavigate () {}

bool
vtkKWScubaToolNavigate::SuspendPickEvents () {

  // Suspend pick events when the button is down so our click-and-drag
  // navigation is fast.
  return this->IsButtonDown();
}

void
vtkKWScubaToolNavigate::DoMouseDown ( vtkKWScubaWindow*,
                                      vtkKWScubaView* iView,
                                      vtkKWScubaLayer*,
                                      float[3] ) {
  iView->StartFastMode();
  iView->GetRenderWindow()->SetDesiredUpdateRate( 30 );
}

void
vtkKWScubaToolNavigate::DoMouseUp ( vtkKWScubaWindow*,
                                    vtkKWScubaView* iView,
                                    vtkKWScubaLayer*,
                                    float[3] ) {
  iView->StopFastMode();
  iView->GetRenderWindow()->SetDesiredUpdateRate( 0 );
}

void
vtkKWScubaToolNavigate::DoMouseDrag ( vtkKWScubaWindow*,
                                      vtkKWScubaView* iView,
                                      vtkKWScubaLayer*,
                                      float[3] ) {

  int delta[2];
  delta[0] = this->GetButtonDeltaX();
  delta[1] = this->GetButtonDeltaY();

  switch( iView->GetDisplayMode() ) {
  case vtkKWScubaView::TwoDee:
    if ( 1 == this->WhichButtonDown() ) {
      iView->PanBetweenWindowCoords( delta );
    } else if ( 2 == this->WhichButtonDown() ) {
      iView->ScrollBetweenWindowCoords( delta );
    } else if ( 3 == this->WhichButtonDown() ) {
      iView->ZoomBetweenWindowCoords( delta );
    }
    break;
  case vtkKWScubaView::ThreeDee:
    if ( 1 == this->WhichButtonDown() ) {
      iView->PanBetweenWindowCoords( delta );
    } else if ( 2 == this->WhichButtonDown() ) {
      iView->RotateBetweenWindowCoords( delta );
    } else if ( 3 == this->WhichButtonDown() ) {
      iView->ZoomBetweenWindowCoords( delta );
    }
    break;
  default:
    break;
  }
}
