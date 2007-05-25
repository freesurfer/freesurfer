/**
 * @file  vtkKWScubaApp.cxx
 * @brief Main app for Scuba
 *
 * The main app code. Parse command line options and does some initial
 * set up.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: kteich $
 *    $Date: 2007/05/25 18:18:05 $
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


#include "vtkKWScubaApp.h"
#include "vtkKWScubaWindow.h"
#include "vtkObjectFactory.h"
#include "vtkKWOptionDataBase.h"
#include "IconLoader.h"

using namespace std;

vtkStandardNewMacro( vtkKWScubaApp );
vtkCxxRevisionMacro( vtkKWScubaApp, "$Revision: 1.3 $" );

vtkKWScubaApp::vtkKWScubaApp () :
    mWindow( NULL ) {

  this->GetOptionDataBase()->
    AddEntry( "vtkKWWidget", "SetBackgroundColor", "1 1 1" );

  // Init the icon loader with the app and load our icons.
  try {
    IconLoader::Initialize( this );

    try {
      IconLoader::LoadIconsFromFile( "./ScubaIcons.txt" );
    }
    catch(...) {
      char* pfnFreesurferDir = getenv( "FREESURFER_HOME" );
      if( NULL != pfnFreesurferDir ) {
	string fnIcons = 
	  string(pfnFreesurferDir) + "/lib/resource/ScubaIcons.txt";
	IconLoader::LoadIconsFromFile( fnIcons.c_str() );
      }
    }
  }
  catch( exception& e ) {
    cerr << "Error loading icons: " << e.what() << endl;
  }

  // Set our app name and load our settings. This needs to be done
  // before we make our window because the window will try to get
  // stuff from the registry too.
  this->SetName( "Scuba" );
  this->RestoreApplicationSettingsFromRegistry();

  // Create the main window.
  mWindow = vtkKWScubaWindow::New();
  mWindow->SetApplication( this );
  this->AddWindow( mWindow );
  mWindow->Create();

}

vtkKWScubaApp::~vtkKWScubaApp () {

  mWindow->Delete();

  this->SaveApplicationSettingsToRegistry();
}

void
vtkKWScubaApp::Start ( int argc, char* argv[] ) {

  mWindow->Display();

  // Command line arguments.
  for ( int i = 1; i < argc; i++ ) {

    // If this isn't a dash-option, continue.
    if ( argv[i] && *argv[i] != '-' ) {
      cerr << "Unrecognized option: " << argv[i] << endl;
    }

    // Now that we're sure it's a dash-option, remove the dashes.
    string arg( argv[i] );
    if ( arg[0] == '-' ) arg = arg.replace( 0, 1, "" );
    if ( arg[0] == '-' ) arg = arg.replace( 0, 1, "" );

    // Check the possibilites.
    if ( arg == "v" || arg == "volume") {

      // Volume. All following arguments are volume names until the
      // next dash. Keep loading volumes until we see another -. For
      // each one, try to load it. If failed, put up an error
      // message and try the next one.
      do {
        string s( argv[++i] );
        try {
          LoadVolume( s.c_str() );
        } catch ( exception& e ) {
          cerr << "Error loading volume \"" << s << "\": " << e.what();
        }
      } while ( i+1 < argc && argv[i+1][0] != '-' );

    }

    if ( arg == "s" || arg == "surface") {

      // Surface. All following arguments are surface names until the
      // next dash. Keep loading surfaces until we see another -. For
      // each one, try to load it. If failed, put up an error
      // message and try the next one.
      do {
        string s( argv[++i] );
        try {
          LoadSurface( s.c_str() );
        } catch ( exception& e ) {
          cerr << "Error loading surface \"" << s << "\": " << e.what();
        }
      } while ( i+1 < argc && argv[i+1][0] != '-' );

    }

    if ( arg == "d" || arg == "dti") {

      // DTI. All following arguments are DTIs names until the
      // next dash. Keep loading DTIs until we see another -. For
      // each one, try to load it. If failed, put up an error
      // message and try the next one.
      do {
        string s( argv[++i] );
        try {
          LoadDTI( s.c_str() );
        } catch ( exception& e ) {
          cerr << "Error loading DTI \"" << s << "\": " << e.what();
        }
      } while ( i+1 < argc && argv[i+1][0] != '-' );

    }
  }

  this->Superclass::Start( argc, argv );
}

void
vtkKWScubaApp::RestoreApplicationSettingsFromRegistry () {

  vtkKWApplication::RestoreApplicationSettingsFromRegistry();
}

void
vtkKWScubaApp::SaveApplicationSettingsToRegistry () {

  vtkKWApplication::SaveApplicationSettingsToRegistry();
}

void
vtkKWScubaApp::LoadVolume ( const char* ifnVolume ) {

  // Just pass to the window.
  if ( mWindow )
    mWindow->LoadVolume( ifnVolume );
}

void
vtkKWScubaApp::LoadSurface ( const char* ifnSurface ) {

  // Just pass to the window.
  if ( mWindow )
    mWindow->LoadSurface( ifnSurface );
}

void
vtkKWScubaApp::LoadDTI ( const char* ifnDTI ) {

  // Just pass to the window.
  if ( mWindow )
    mWindow->LoadDTI( ifnDTI );
}

void
vtkKWScubaApp::LoadPath ( const char* ifnPath ) {

  // Just pass to the window.
  if ( mWindow )
    mWindow->LoadPath( ifnPath );
}
