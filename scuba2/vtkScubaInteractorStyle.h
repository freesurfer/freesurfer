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
