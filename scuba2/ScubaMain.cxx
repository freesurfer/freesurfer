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
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:40 $
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


#include <string>
#include <vtksys/SystemTools.hxx>
#include <vtksys/CommandLineArguments.hxx>

#include "vtkKWScubaApp.h"

extern "C" {

  extern int Scubalib_SafeInit( Tcl_Interp* );

}

using namespace std;

const char* Progname = "scuba2";

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
