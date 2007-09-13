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
 *    $Author: kteich $
 *    $Date: 2007/09/13 20:58:21 $
 *    $Revision: 1.4 $
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
