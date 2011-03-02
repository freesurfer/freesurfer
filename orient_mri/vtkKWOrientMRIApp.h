/**
 * @file  vtkKWOrientMRIApp.h
 * @brief Command line parsing and startup
 *
 * Creates our window, does some setup stuff, parses command line
 * args.
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


#ifndef __vtkKWOrientMRIApp_h
#define __vtkKWOrientMRIApp_h

#include <string>

#include "vtkKWApplication.h"
#include "vtkSmartPointer.h"

class vtkKWOrientMRIWindow;

class vtkKWOrientMRIApp : public vtkKWApplication {

public:

  static vtkKWOrientMRIApp* New ();
  vtkTypeRevisionMacro ( vtkKWOrientMRIApp, vtkKWApplication );

  // Description:
  // Override to show our window.
  virtual void Start ( int argc, char* argv[] );

  // Description:
  // Load a volume and display it.
  void LoadVolume ( const char* ifnVolume );

  // Description:
  // Override error message routine to display a dialog box as well.
  virtual void ErrorMessage ( const char* isMessage );

protected:

  vtkKWOrientMRIApp ();
  ~vtkKWOrientMRIApp ();

  //BTX
  // Our window object.
  vtkSmartPointer<vtkKWOrientMRIWindow> mWindow;
  //ETX
};

#endif
