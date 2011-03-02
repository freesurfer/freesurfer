/**
 * @file  vtkKWOrientMRIApp.cxx
 * @brief Command line parsing and startup
 *
 * Creates our window, does some setup stuff, parses command line
 * args.
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

#include <assert.h>

#include "vtkKWOrientMRIApp.h"

#include "IconLoader.h"
#include "vtkKWMessageDialog.h"
#include "vtkKWOptionDataBase.h"
#include "vtkKWOrientMRIWindow.h"
#include "vtkObjectFactory.h"
#include "vtksys/CommandLineArguments.hxx"
#include "vtksys/SystemTools.hxx"

using namespace std;

vtkStandardNewMacro( vtkKWOrientMRIApp );
vtkCxxRevisionMacro( vtkKWOrientMRIApp, "$Revision: 1.7 $" );

vtkKWOrientMRIApp::vtkKWOrientMRIApp () {

  this->GetOptionDataBase()->
    AddEntry( "vtkKWWidget", "SetBackgroundColor", "1 1 1" );

  // Init the icon loader with the app and load our icons.
  try {
    IconLoader::Initialize( this );

    try {
      IconLoader::LoadIconsFromFile( "./OrientMRIIcons.txt" );
    }
    catch(...) {
      char* pfnFreesurferDir = getenv( "FREESURFER_HOME" );
      if( NULL != pfnFreesurferDir ) {
	string fnIcons = 
	  string(pfnFreesurferDir) + "/lib/resource/OrientMRIIcons.txt";
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
  this->SetName( "orient_mri" );
  this->RestoreApplicationSettingsFromRegistry();

  // Create the main window.
  mWindow = vtkSmartPointer<vtkKWOrientMRIWindow>::New();
  mWindow->SetApplication( this );
  this->AddWindow( mWindow );
  mWindow->Create();
}

vtkKWOrientMRIApp::~vtkKWOrientMRIApp () {

}

void
vtkKWOrientMRIApp::Start ( int argc, char* argv[] ) {

  assert( mWindow.GetPointer() );

  this->SetHelpDialogStartingPage("https://surfer.nmr.mgh.harvard.edu/fswiki/orient_5fmri");

  mWindow->Display();

  // Create our command line argument parser.
  vtksys::CommandLineArguments args;
  args.Initialize( argc, argv );

  // Add the arguments we'll look for.
  string fnVolume = "";
  args.AddArgument( "--volume", args.SPACE_ARGUMENT, &fnVolume,
		    "A volume file to load" );

  // Try and parse the arguments. If there was an error, print our
  // help message and quit.
  if( !args.Parse() ) {
    cerr << "Error parsing arguments." << endl;
    cerr << args.GetHelp() << endl;
    exit( 1 );
  }

  // If they gave us a volume name, load it.
  if( fnVolume != "" )
    this->LoadVolume( fnVolume.c_str() );


  this->Superclass::Start( argc, argv );
}

void
vtkKWOrientMRIApp::LoadVolume ( const char* ifnVolume ) {

  assert( mWindow.GetPointer() );

  mWindow->LoadVolume( ifnVolume );
}

void
vtkKWOrientMRIApp::ErrorMessage ( const char* isMessage ) {

  vtkSmartPointer<vtkKWMessageDialog> dialog = 
    vtkSmartPointer<vtkKWMessageDialog>::New();
  dialog->SetStyleToMessage();
  dialog->SetOptions( vtkKWMessageDialog::ErrorIcon );
  dialog->SetApplication( this );
  dialog->Create();
  dialog->SetText( isMessage );
  dialog->Invoke();

  this->vtkKWApplication::ErrorMessage( isMessage );
}

