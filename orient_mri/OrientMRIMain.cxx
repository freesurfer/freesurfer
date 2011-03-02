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
 *    $Date: 2011/03/02 00:04:35 $
 *    $Revision: 1.7 $
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
