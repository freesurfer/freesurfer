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
  if( !interp ) {
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
