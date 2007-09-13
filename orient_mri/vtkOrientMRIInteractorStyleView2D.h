/**
 * @file  vtkOrientMRIInteractorStyleView2D.h
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

#ifndef vtkOrientMRIInteractorStyleView2D_h
#define vtkOrientMRIInteractorStyleView2D_h

#include "vtkInteractorStyleTrackballCamera.h"
#include "vtkSmartPointer.h"

class vtkWindowLevelLookupTable;

class vtkOrientMRIInteractorStyleView2D : public vtkInteractorStyleTrackballCamera {

public:
  static vtkOrientMRIInteractorStyleView2D *New();
  vtkTypeRevisionMacro(vtkOrientMRIInteractorStyleView2D,vtkInteractorStyle);

  // Description:
  // Sets the lookup table to do window/level on.
  void SetWindowLevelTable ( vtkWindowLevelLookupTable* iTable );

  // Description:
  // Override left button stuff to do the window/level operation.
  virtual void OnMouseMove ();
  virtual void OnLeftButtonDown ();
  virtual void OnLeftButtonUp ();

  // Description:
  // Do a window/level operation on the lookup table we have.
  void WindowLevel ();

protected:
  vtkOrientMRIInteractorStyleView2D();
  ~vtkOrientMRIInteractorStyleView2D();

  // Start and end our window/level operation.
  void StartWindowLevel ();
  void EndWindowLevel ();

  // Our custom state.
  enum { IS_WindowLevel = 1000 };

  // The window/level color table we'll modify with the left button.
  vtkSmartPointer<vtkWindowLevelLookupTable> LookupTable;

  // Save original window info.
  double OriginalWindow;

private:
  vtkOrientMRIInteractorStyleView2D(const vtkOrientMRIInteractorStyleView2D&);
  void operator=(const vtkOrientMRIInteractorStyleView2D&);
};

#endif
