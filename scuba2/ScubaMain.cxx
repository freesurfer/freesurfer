/**
 * @file  ScubaMain.cxx
 * @brief Main file for Scuba
 *
 * Simple main file that inits the library and starts and runs the
 * application.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: kteich $
 *    $Date: 2007/04/06 22:23:04 $
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


#include <string>
#include <vtksys/SystemTools.hxx>
#include <vtksys/CommandLineArguments.hxx>

#include "vtkKWScubaApp.h"

extern "C" {

  extern int Scubalib_SafeInit( Tcl_Interp* );

}

using namespace std;

char* Progname = "scuba2";

int main ( int argc, char** argv ) {

  // Initialize Tcl.
  Tcl_Interp* interp = vtkKWApplication::InitializeTcl( argc, argv, &cerr );
  if ( !interp ) {
    cerr << "Error initializing Tcl." << endl;
    return 1;
  }

  // Init our Tcl wrapping code.
  Scubalib_SafeInit( interp );

  // Create and run the app.
  vtkKWScubaApp* app = vtkKWScubaApp::New();
  app->Start( argc, argv );
  int rApp = app->GetExitStatus();

  // Delete.
  app->Delete();

  return rApp;
}
