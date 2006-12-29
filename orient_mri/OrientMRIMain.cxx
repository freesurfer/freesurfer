/**
 * @file  OrientMRIMain.cxx
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 02:09:11 $
 *    $Revision: 1.3 $
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
#include "OrientMRICustomIcons.h"

extern "C" {

  extern int Orientmrilib_SafeInit( Tcl_Interp* );

}

using namespace std;

char* Progname = "orient_mri";

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
  vtkKWOrientMRIApp* app = vtkKWOrientMRIApp::New();

  // Init the custom icons lib.
  OrientMRICustomIcons::Initialize( app );

  // Run the app.
  app->Start( argc, argv );
  int rApp = app->GetExitStatus();

  // Delete.
  app->Delete();

  OrientMRICustomIcons::ShutDown();

  return rApp;
}
