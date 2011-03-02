/**
 * @file  vtkScubaInteractorStyle.h
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


// .NAME vtkScubaInteractorStyle - Pass events to vtkKWScubaView
// .SECTION Description
// This is an implementation of the VTK interactor style method, a
// listener that responds to interaction events. This class takes
// those events and passes them to the vtkKWScubaView.

#ifndef vtkScubaInteractorStyle_h
#define vtkScubaInteractorStyle_h

#include "vtkInteractorStyle.h"

class vtkKWScubaWindow;

class vtkScubaInteractorStyle : public vtkInteractorStyle {

public:

  static vtkScubaInteractorStyle* New();
  vtkTypeRevisionMacro( vtkScubaInteractorStyle, vtkInteractorStyle );

  // Description:
  // Set the window. When an event happens, the window will be passed
  // along to the view with the event call. It will also tell the
  // window when an event is done being handled.
  void SetWindow ( vtkKWScubaWindow* iWindow );

  // Description:
  // Implement vtkInteractorStyle callbacks so we can pass events our
  // way.
  virtual void OnMouseMove ();
  virtual void OnLeftButtonDown ();
  virtual void OnLeftButtonUp ();
  virtual void OnMiddleButtonDown ();
  virtual void OnMiddleButtonUp ();
  virtual void OnRightButtonDown ();
  virtual void OnRightButtonUp ();
  virtual void OnKeyDown ();
  virtual void OnKeyUp ();
  virtual void OnEnter ();
  virtual void OnLeave ();

protected:

  vtkScubaInteractorStyle();
  virtual ~vtkScubaInteractorStyle();

  vtkKWScubaWindow* mWindow;

private:
  vtkScubaInteractorStyle( const vtkScubaInteractorStyle& );
  void operator=( const vtkScubaInteractorStyle& );
};


#endif
