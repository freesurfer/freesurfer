/**
 * @file  IconLoaderTest.cxx
 * @brief Tests IconLoader class
 *
 */
/*
 * Original Author: Nick Schmansky
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2008/06/17 21:34:27 $
 *    $Revision: 1.2 $
 *
 * Copyright (C) 2008,
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

#include <iostream>
#include <stdexcept>

#include "vtkObjectFactory.h"
#include "vtkKWApplication.h"
#include "IconLoader.h"
#include "IconLoaderTest.h"

extern "C" {

#include "unistd.h" // getcwd
#include "tix.h"
  
  extern int Qdeclib_SafeInit( Tcl_Interp* iInterp );
  extern int Blt_Init( Tcl_Interp* iInterp );

}

using namespace std;

vtkStandardNewMacro( IconLoaderTest );

IconLoaderTest::IconLoaderTest () :
  vtkKWApplication() {

  // Init the icon loader with the app and load our icons.
  try {
    IconLoader::Initialize( this );
    
    try {
      IconLoader::LoadIconsFromFile( "./IconLoaderTestIcons.txt" );
    }
    catch(...) {
      char* pfnFreesurferDir = getenv( "FREESURFER_HOME" );
      if( NULL != pfnFreesurferDir ) {
        string fnIcons =
          string(pfnFreesurferDir) + "/lib/resource/QdecIcons.txt";
        IconLoader::LoadIconsFromFile( fnIcons.c_str() );
      }
    }
  }
  catch( exception& e ) {
    cerr << "Error loading icons: " << e.what() << endl;
  }
  cout << "Success loading icons" << endl;
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

  return 0;
}
