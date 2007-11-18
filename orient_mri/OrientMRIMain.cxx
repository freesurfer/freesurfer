/**
 * @file  OrientMRIMain.cxx
 * @brief main function
 *
 * Starts up our libraries and runs the application.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2007/11/18 03:03:37 $
 *    $Revision: 1.6 $
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


#include <string>
#include <vtksys/SystemTools.hxx>
#include <vtksys/CommandLineArguments.hxx>

#include "vtkKWOrientMRIApp.h"
#include "vtkSmartPointer.h"

extern "C" {

  extern int Orientmrilib_SafeInit( Tcl_Interp* );

}

using namespace std;

const char* Progname = "orient_mri";

int main ( int argc, char** argv ) {

  // Initialize Tcl.
  Tcl_Interp* interp = vtkKWApplication::InitializeTcl( argc, argv, &cerr );
  if ( !interp ) {
    cerr << "Error initializing Tcl." << endl;
    return 1;
  }

  // Init our Tcl wrapping code.
  Orientmrilib_SafeInit( interp );

  // Create the app.
  vtkSmartPointer<vtkKWOrientMRIApp> app =
    vtkSmartPointer<vtkKWOrientMRIApp>::New();

  // Run the app.
  app->Start( argc, argv );
  int rApp = app->GetExitStatus();

  return rApp;
}
