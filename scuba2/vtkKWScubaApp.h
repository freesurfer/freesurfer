/**
 * @file  vtkKWScubaApp.h
 * @brief Main app for Scuba
 *
 * The main app code. Parse command line options and does some initial
 * set up.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:40 $
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


// .NAME vtkKWScubaApp - Scuba application
// .SECTION Description
// Starts up the up, creates a vtkKWScubaWindow, and does some command
// line option processing.

#ifndef __vtkKWScubaApp_h
#define __vtkKWScubaApp_h

#include <string>
#include "vtkKWApplication.h"
#include "vtkSmartPointer.h"

class vtkKWScubaWindow;

class vtkKWScubaApp : public vtkKWApplication {

public:

  static vtkKWScubaApp* New ();
  vtkTypeRevisionMacro ( vtkKWScubaApp, vtkKWApplication );

  // Description:
  // Override to show our window.
  virtual void Start ( int argc, char* argv[] );

  // Description:
  // Tells the window to load something.
  void LoadVolume ( const char* ifnVolume );
  void LoadSurface ( const char* ifnSurface );
  void LoadDTI ( const char* ifnDTI );
  void LoadPath ( const char* ifnPath );

protected:

  vtkKWScubaApp ();
  ~vtkKWScubaApp ();

  //BTX

  // Description:
  // If we have a subjects dir and subject name, and this is not a
  // full path file name, create a full path using the subjects dir
  // and the given subdir.
  std::string FormatFileNameUsingSubjectName ( const char* ifnMain, 
					       const char* ifnSubDir );

  // Pointer to our main window.
  vtkSmartPointer<vtkKWScubaWindow> mWindow;

  // Our subject name, if set.
  std::string msSubjectName;
  //ETX
};

#endif
