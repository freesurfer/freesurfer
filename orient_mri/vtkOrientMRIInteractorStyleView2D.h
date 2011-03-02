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
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:35 $
 *    $Revision: 1.3 $
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
  // Set the ortho points we're working on. 
  void SetPointsAndLines ( double* iPointsButton1Start, 
			   double* iPointsButton1End,
			   double* iPointsButton2Start, 
			   double* iPointsButton2End,
			   double* iPointsButton3Start, 
			   double* iPointsButton3End );

  // Description:
  // Override left button stuff to do the window/level
  // operation. Override control-clicks on all buttons to do ortho
  // line drawing.
  virtual void OnMouseMove ();
  virtual void OnLeftButtonDown ();
  virtual void OnLeftButtonUp ();
  virtual void OnMiddleButtonDown ();
  virtual void OnMiddleButtonUp ();
  virtual void OnRightButtonDown ();
  virtual void OnRightButtonUp ();

  // Description:
  // Do a window/level operation on the lookup table we have.
  void WindowLevel ();

  // Description:
  // Do a ortho line edit operation on the points and lines of the
  // index ActiveOrthoLine.
  void EditOrthoLine ();

protected:
  vtkOrientMRIInteractorStyleView2D();
  ~vtkOrientMRIInteractorStyleView2D();

  // Start and end our window/level operation.
  void StartWindowLevel ();
  void EndWindowLevel ();

  // Start and end our ortho line editing operation.
  void StartOrthoLineEdit ( int iWhich );
  void EndOrthoLineEdit ();

  // Our custom states.
  enum { IS_WindowLevel = 1000, IS_OrthoLineEdit };

  // The window/level color table we'll modify with the left button.
  vtkSmartPointer<vtkWindowLevelLookupTable> LookupTable;

  // The points we'll work with.
  double* Points[3][2];

  // Save original window info.
  double OriginalWindow;

  // Which ortho line we're drawing.
  int ActiveOrthoLine;

  // The first position of the line start (gotten at mouse down).
  double ActiveOrthoLineStart[3];

private:
  vtkOrientMRIInteractorStyleView2D(const vtkOrientMRIInteractorStyleView2D&);
  void operator=(const vtkOrientMRIInteractorStyleView2D&);
};

#endif
