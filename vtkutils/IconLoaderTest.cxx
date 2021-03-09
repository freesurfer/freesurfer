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

#include "vtkObjectFactory.h"
#include "vtkKWApplication.h"
#include "IconLoader.h"
#include "IconLoaderTest.h"

extern "C" {
#include "unistd.h" // getcwd
#include "tix.h"
    extern int Blt_Init( Tcl_Interp* iInterp );
}

using namespace std;

static int errs=0;

vtkStandardNewMacro( IconLoaderTest );

IconLoaderTest::IconLoaderTest () :
  vtkKWApplication() {

  // Init the icon loader with the app and load our icons.
  try {
    IconLoader::Initialize( this );
    
    try {
      errs += IconLoader::LoadIconsFromFile( "./IconLoaderTestIcons.txt" );
    }
    catch(...) {
      char* pfnFreesurferDir = getenv( "FREESURFER_HOME" );
      if( NULL != pfnFreesurferDir ) {
        string fnIcons =
          string(pfnFreesurferDir) + "/lib/resource/QdecIcons.txt";
        errs += IconLoader::LoadIconsFromFile( fnIcons.c_str() );
      }
    }
  }
  catch( exception& e ) {
    cerr << "Error loading icons: " << e.what() << endl;
    errs++;
  }

  if (0 == errs) cout << "Success loading icons" << endl;
};

IconLoaderTest::~IconLoaderTest () {

  IconLoader::ShutDown();
}

int main ( int argc, char** argv ) {

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

  IconLoaderTest* app = IconLoaderTest::New();
  app->Delete();

  return errs;
}
