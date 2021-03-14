/**
 * @brief Tests IconLoader class
 *
 */
/*
 * Original Author: Nick Schmansky
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

#include <iostream>
#include <stdexcept>

#include "IconLoader.h"
#include "IconLoaderTest.h"
#include "vtkKWApplication.h"
#include "vtkObjectFactory.h"

extern "C" {
#include "tix.h"
#include "unistd.h" // getcwd
extern int Blt_Init(Tcl_Interp *iInterp);
}

static int errs = 0;

vtkStandardNewMacro(IconLoaderTest);

IconLoaderTest::IconLoaderTest() : vtkKWApplication() {

  // Init the icon loader with the app and load our icons.
  try {
    IconLoader::Initialize(this);

    try {
      errs += IconLoader::LoadIconsFromFile("./IconLoaderTestIcons.txt");
    } catch (...) {
      char *pfnFreesurferDir = getenv("FREESURFER_HOME");
      if (NULL != pfnFreesurferDir) {
        std::string fnIcons =
            string(pfnFreesurferDir) + "/lib/resource/QdecIcons.txt";
        errs += IconLoader::LoadIconsFromFile(fnIcons.c_str());
      }
    }
  } catch (exception &e) {
    std::cerr << "Error loading icons: " << e.what() << std::endl;
    errs++;
  }

  if (0 == errs)
    std::cout << "Success loading icons" << std::endl;
};

IconLoaderTest::~IconLoaderTest() { IconLoader::ShutDown(); }

int main(int argc, char **argv) {

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

  IconLoaderTest *app = IconLoaderTest::New();
  app->Delete();

  return errs;
}
