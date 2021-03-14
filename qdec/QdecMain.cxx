/**
 * @brief Starts the application
 *
 * This starts the Tcl interpreter, inits our application library, and
 * starts up the application object.
 */
/*
 * Original Author: Kevin Teich
 *
 * Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
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
#include <vtksys/CommandLineArguments.hxx>
#include <vtksys/SystemTools.hxx>

#include "vtkKWQdecApp.h"

#include "diag.h"
#include "error.h"
#include "tix.h"
#include "unistd.h" // getcwd

extern int Qdeclib_SafeInit(Tcl_Interp *iInterp);
extern int Blt_Init(Tcl_Interp *iInterp);

const char *Progname = "qdec";

int main(int argc, char **argv) {

  ErrorInit(NULL, NULL, NULL);
  DiagInit(NULL, NULL, NULL);

  // if SUBJECTS_DIR is not set, then set it to the current working dir
  if (NULL == getenv("SUBJECTS_DIR")) {
    if (setenv("SUBJECTS_DIR", getcwd(NULL, 0), 1)) {
      std::cerr << std::endl
                << "ERROR: failure setting SUBJECTS_DIR to cwd." << std::endl;
      return 1;
    }
  }

  // Initialize Tcl.
  Tcl_Interp *interp = vtkKWApplication::InitializeTcl(argc, argv, &cerr);
  if (!interp) {
    std::cerr << "Error initializing Tcl." << std::endl;
    return 1;
  }

  // Init Tix manually.
  int rTcl = Tix_Init(interp);
  if (TCL_OK != rTcl) {
    const char *sResult = Tcl_GetStringResult(interp);
    std::cerr << "Tix_Init returned not TCL_OK: " << sResult << std::endl;
    return 1;
  }

  // Init Blt manually.
  rTcl = Blt_Init(interp);
  if (TCL_OK != rTcl) {
    const char *sResult = Tcl_GetStringResult(interp);
    std::cerr << "Blt_Init returned not TCL_OK: " << sResult << std::endl;
    return 1;
  }

  // Init our Tcl wrapping code.
  Qdeclib_SafeInit(interp);

  // Create the app.
  vtkKWQdecApp *app = vtkKWQdecApp::New();

  // Run the app.
  app->Start(argc, argv);
  int rApp = app->GetExitStatus();

  // Delete.
  app->Delete();

  return rApp;
}
