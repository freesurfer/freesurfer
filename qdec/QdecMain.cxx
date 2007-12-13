/**
 * @file  QdecMain.cxx
 * @brief Starts the application
 *
 * This starts the Tcl interpreter, inits our application library, and
 * starts up the application object.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2007/12/13 22:22:16 $
 *    $Revision: 1.1 $
 *
 * Copyright (C) 2007,
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

#include "vtkKWQdecApp.h"

extern "C" {

#include "unistd.h" // getcwd
#include "tix.h"
  
  extern int Qdeclib_SafeInit( Tcl_Interp* iInterp );
  extern int Blt_Init( Tcl_Interp* iInterp );

}

using namespace std;

const char* Progname = "qdec";

int main ( int argc, char** argv ) {

  // if SUBJECTS_DIR is not set, then set it to the current working dir
  if ( NULL == getenv( "SUBJECTS_DIR" ) ) {
    if ( setenv( "SUBJECTS_DIR", getcwd(NULL,0), 1) ) {
      cerr << endl << "ERROR: failure setting SUBJECTS_DIR to cwd." << endl;
      return 1;
    }
  }

  // Initialize Tcl.
  Tcl_Interp* interp = vtkKWApplication::InitializeTcl( argc, argv, &cerr );
  if ( !interp ) {
    cerr << "Error initializing Tcl." << endl;
    return 1;
  }

  // Init Tix manually.
  int rTcl = Tix_Init( interp );
  if ( TCL_OK != rTcl ) {
    const char* sResult = Tcl_GetStringResult( interp );
    cerr <<  "Tix_Init returned not TCL_OK: " << sResult << endl;
    return 1;
  }

  // Init Blt manually.
  rTcl = Blt_Init( interp );
  if ( TCL_OK != rTcl ) {
    const char* sResult = Tcl_GetStringResult( interp );
    cerr <<  "Blt_Init returned not TCL_OK: " << sResult << endl;
    return 1;
  }

  // Init our Tcl wrapping code.
  Qdeclib_SafeInit( interp );

  // Create the app.
  vtkKWQdecApp* app = vtkKWQdecApp::New();

  // Run the app.
  app->Start( argc, argv );
  int rApp = app->GetExitStatus();

  // Delete.
  app->Delete();

  return rApp;
}
