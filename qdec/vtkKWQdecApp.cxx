/**
 * @brief Command line parsing and registry
 *
 * Application code that parses the command line options and passes
 * them to the main window. Also gets registry options pertaining to
 * the window size.
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

#define _GLIBCXX_USE_CXX11_ABI 0

#include <stdexcept>

#include "IconLoader.h"
#include "vtkKWQdecApp.h"
#include "vtkKWQdecWindow.h"
#include "vtkObjectFactory.h"
#include "vtkKWDialog.h"
#include "vtkKWMessageDialog.h"
#include "vtkKWOptionDataBase.h"
#include "vtkKWPushButton.h"
#include "vtkKWSplashScreen.h"
#include "vtkKWText.h"
#include "vtkKWTextWithScrollbars.h"
#include "vtkSmartPointer.h"
#include "vtksys/CommandLineArguments.hxx"
#include "vtksys/SystemTools.hxx"
#include "QdecUtilities.h"
#include "FsgdfPlot.h"

using namespace std;

vtkStandardNewMacro( vtkKWQdecApp );

vtkKWQdecApp::vtkKWQdecApp () :
  vtkKWApplication() {

  // Set us up with some global styles.
  this->GetOptionDataBase()->
    AddEntry( "vtkKWWidget", "SetBackgroundColor", "1 1 1" );
  this->GetOptionDataBase()->
    AddFontOptions( "helvetica 12" );

  // Init the icon loader with the app and load our icons.
  try {
    IconLoader::Initialize( this );

    try {
      IconLoader::LoadIconsFromFile( "./QdecIcons.txt" );
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

  // Set some application stuff.
  this->SetName( "Qdec" );
  this->SetMajorVersion( 1 ); // v1.5
  this->SetMinorVersion( 5 );
  this->SetHelpDialogStartingPage
    ("https://surfer.nmr.mgh.harvard.edu/fswiki/Qdec");
  this->SupportSplashScreenOn ( );
  this->SplashScreenVisibilityOn ( );

  // Set the window size from the registry values.
  this->RestoreApplicationSettingsFromRegistry();

}

vtkKWQdecApp::~vtkKWQdecApp () {

  this->SaveApplicationSettingsToRegistry();

  IconLoader::ShutDown();
}

void
vtkKWQdecApp::Start ( int argc, char* argv[] ) {

  vtkKWSplashScreen* pSplash = this->GetSplashScreen();
  char* pfnFreesurferDir = getenv( "FREESURFER_HOME" );
  if( NULL != pfnFreesurferDir ) {
    string fSplash = string(pfnFreesurferDir) + "/lib/images/qdec_splash.png";
    pSplash->ReadImage( fSplash.c_str() );
    pSplash->Display();
  }

  // Create our command line argument parser.
  vtksys::CommandLineArguments args;
  args.Initialize( argc, argv );

  // Add the arguments we'll look for.
  bool bHelp = false;
  args.AddArgument( "--help", args.NO_ARGUMENT, &bHelp,
                    "Print option information" );

  string fnTable;
  args.AddArgument( "--table", args.SPACE_ARGUMENT, &fnTable,
                    "Data table to load" );

  string fnProject;
  args.AddArgument( "--project", args.SPACE_ARGUMENT, &fnProject,
                    "Project file (*.qdec) to load" );

  vector<string> lfnSurfaces;
#if ((VTK_MAJOR_VERSION == 5)&&(VTK_MINOR_VERSION > 0))
  args.AddArgument( "--surface", args.MULTI_ARGUMENT, &lfnSurfaces,
                    "A surface file or list of files to load" );
#else
  string fnSurface;
  args.AddArgument( "--surface", args.SPACE_ARGUMENT, &fnSurface,
                    "A surface file to load" );
#endif

  string fnGDF;
  args.AddArgument( "--gdf", args.SPACE_ARGUMENT, &fnGDF,
                    "FSGDF file to load" );
  args.AddArgument( "--fsgd", args.SPACE_ARGUMENT, &fnGDF,
                    "FSGDF file to load" );

  vector<string> lfnScalars;
#if ((VTK_MAJOR_VERSION == 5)&&(VTK_MINOR_VERSION > 0))
  args.AddArgument( "--scalar", args.MULTI_ARGUMENT, &lfnScalars,
                    "A scalar file or list of files to load" );
#else
  string fnScalar;
  args.AddArgument( "--scalar", args.SPACE_ARGUMENT, &fnScalar,
                    "A scalar file to load" );
#endif

  string fnCurvature;
  args.AddArgument( "--curvature", args.SPACE_ARGUMENT, &fnCurvature,
                    "Curvature file to load" );

  string fnAnnotation;
  args.AddArgument( "--annotation", args.SPACE_ARGUMENT, &fnAnnotation,
                    "Annotation file to load" );

  vector<string> lfnOverlay;
#if ((VTK_MAJOR_VERSION == 5)&&(VTK_MINOR_VERSION > 0))
  args.AddArgument( "--overlay", args.MULTI_ARGUMENT, &lfnOverlay,
                    "Overlay file and color file to load" );
#endif

  string fnLabel;
  args.AddArgument( "--label", args.SPACE_ARGUMENT, &fnLabel,
                    "Label file to load" );

  string sAverageSubject;
  args.AddArgument( "--average", args.SPACE_ARGUMENT, &sAverageSubject,
                    "Name of average (common-space) subject" );

  // Try and parse the arguments. If there was an error, print our
  // help message and quit.
  if( !args.Parse() ) {
    cerr << "Error parsing arguments." << endl;
    cerr << args.GetHelp() << endl;
    exit( 1 );
  }
  if( bHelp ) {
    cerr << "USAGE: qdec [options...]" << endl << endl;
    cerr << "Options:" << endl;
    cerr << args.GetHelp() << endl;
    exit( 1 );
  }

#if ((VTK_MAJOR_VERSION == 5)&&(VTK_MINOR_VERSION > 0))
#else
  if (fnSurface.size()) lfnSurfaces.push_back(fnSurface);
  if (fnScalar.size()) lfnScalars.push_back(fnScalar);
#endif

  // Check results.
  if( lfnOverlay.size() != 0 &&
      lfnOverlay.size() != 2 ) {
    cerr << "Error parsing --overlay option; needs two arguments" << endl;
    cerr << args.GetHelp() << endl;
    exit( 1 );
  }

  // Create the main window.
  mWindow = vtkSmartPointer<vtkKWQdecWindow>::New();
  mWindow->SetApplication( this );
  this->AddWindow( mWindow );

  // Find out if we should use the histogram editor.
  char* pDontUseHistogramEditor = getenv( "QDEC_DONT_USE_HISTOGRAM_EDITOR" );
  if( NULL == pDontUseHistogramEditor )
    mWindow->SetUseHistogramEditor( true );
  else
    mWindow->SetUseHistogramEditor( false );

  // Create the window.
  mWindow->Create();
  mWindow->FinishCreating();

  // Show the window.
  mWindow->Display();

  // Draw the window and stuff before we start processing command line
  // options.
  this->ProcessPendingEvents();

  // Load up the data we got.
  try {

    if( !fnTable.empty() )
      this->LoadDataTable( fnTable.c_str() );

    if( !fnProject.empty() )
      this->LoadProjectFile( fnProject.c_str() );

    vector<string>::iterator tfn;
    for( tfn = lfnSurfaces.begin(); tfn != lfnSurfaces.end(); ++tfn )
      this->LoadSurface( tfn->c_str() );

    if( !fnGDF.empty() )
      this->LoadGDFFile( fnGDF.c_str() );

    for( tfn = lfnScalars.begin(); tfn != lfnScalars.end(); ++tfn )
      this->LoadSurfaceScalars( tfn->c_str() );

    if( !fnCurvature.empty() )
      this->LoadSurfaceCurvatureScalars( fnCurvature.c_str() );

    if( !fnAnnotation.empty() )
      this->LoadAnnotation( fnAnnotation.c_str() );

    if( lfnOverlay.size() == 2 )
      this->LoadSurfaceOverlayScalars( lfnOverlay[0].c_str(),
                                       lfnOverlay[1].c_str() );

    if( !fnLabel.empty() )
      this->LoadLabel( fnLabel.c_str() );

    if( !sAverageSubject.empty() )
      this->SetAverageSubject( sAverageSubject.c_str() );

  }
  catch( exception& e ) {
    this->ErrorMessage( e.what() );
  }

  this->Superclass::Start( argc, argv );
}

int
vtkKWQdecApp::Exit () {

  return this->Superclass::Exit();
}

void
vtkKWQdecApp::AddAboutText( ostream &os) {

  this->Superclass::AddAboutText( os );

  string buildStamp = "Build: ";
  buildStamp += __DATE__ ;
  buildStamp += " " ;
  buildStamp += __TIME__ ;
  buildStamp += "\n  - Copyright (c) 2007-2013\n";
  buildStamp += "    The General Hospital Corporation (Boston, MA),\n";
  buildStamp += "    Martinos Center for Biomedical Imaging,\n";
  buildStamp += "    http://www.nmr.mgh.harvard.edu\n";

  os.write( buildStamp.c_str(), buildStamp.size() );
}

void
vtkKWQdecApp::LoadDataTable ( const char* ifnDataTable ) {
  if ( mWindow.GetPointer() )
    mWindow->LoadDataTable( ifnDataTable );
}

void
vtkKWQdecApp::LoadProjectFile ( const char* ifnProject ) {
  if ( mWindow.GetPointer() )
    mWindow->LoadProjectFile( ifnProject );
}
void
vtkKWQdecApp::LoadSurface ( const char* ifnSurface ) {
  if ( mWindow.GetPointer() ) {

    string fnSurface = ifnSurface;
    string::size_type lastSlash = fnSurface.rfind( "/" );
    string sLabel = fnSurface.substr( lastSlash+1, string::npos ).c_str();

    mWindow->LoadSurface( ifnSurface, sLabel.c_str() );
  }
}

void
vtkKWQdecApp::LoadGDFFile ( const char* ifnGDFFile ) {
  if ( mWindow.GetPointer() )
    mWindow->LoadGDFFile( ifnGDFFile );
}

void
vtkKWQdecApp::LoadSurfaceScalars ( const char* ifnScalars ) {
  if ( mWindow.GetPointer() )
    mWindow->LoadSurfaceScalars( ifnScalars );
}

void
vtkKWQdecApp::LoadSurfaceCurvatureScalars ( const char* ifnScalars ) {
  if ( mWindow.GetPointer() )
    mWindow->LoadSurfaceCurvatureScalars( ifnScalars );
}

void
vtkKWQdecApp::LoadAnnotation ( const char* ifnAnnotation ) {
  if ( mWindow.GetPointer() )
    mWindow->LoadAnnotation( ifnAnnotation );
}

void
vtkKWQdecApp::LoadSurfaceOverlayScalars ( const char* ifnScalars,
                                          const char* ifnColors ) {
  if ( mWindow.GetPointer() )
    mWindow->LoadSurfaceOverlayScalars( ifnScalars, ifnColors );
}

void
vtkKWQdecApp::LoadLabel ( const char* ifnLabel ) {
  if ( mWindow.GetPointer() )
    mWindow->LoadLabel( ifnLabel );
}

void
vtkKWQdecApp::SetAverageSubject ( const char* isAvgSubj ) {
  if ( mWindow.GetPointer() )
    mWindow->SetAverageSubject( isAvgSubj );
}

void
vtkKWQdecApp::ErrorMessage ( const char* isMessage ) {

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

void
vtkKWQdecApp::DisplayHelpDialog ( vtkKWTopLevel* iTop ) {

  if( NULL == mDlogHelp.GetPointer() ) {

    mDlogHelp = vtkSmartPointer<vtkKWDialog>::New();
    mDlogHelp->SetMasterWindow( iTop );
    mDlogHelp->SetApplication( this );
    mDlogHelp->Create();
    mDlogHelp->ModalOff();
    mDlogHelp->SetTitle( "Qdec Help" );
    mDlogHelp->SetSize( 400, 300 );
    mDlogHelp->SetDeleteWindowProtocolCommand( mDlogHelp, "Withdraw" );

    vtkSmartPointer<vtkKWTextWithScrollbars> scrolledText =
      vtkSmartPointer<vtkKWTextWithScrollbars>::New();
    scrolledText->SetParent( mDlogHelp );
    scrolledText->HorizontalScrollbarVisibilityOff();
    scrolledText->Create();

    vtkSmartPointer<vtkKWPushButton> btnClose =
      vtkSmartPointer<vtkKWPushButton>::New();
    btnClose->SetParent( mDlogHelp );
    btnClose->Create();
    btnClose->SetText( "Close" );
    btnClose->SetCommand( mDlogHelp, "Withdraw" );
    btnClose->SetWidth( 10 );

    this->Script( "pack %s -side top -fill both -expand yes -padx 10 -pady 10",
                  scrolledText->GetWidgetName() );
    this->Script( "pack %s -side top -anchor e -padx 10 -pady 10",
                  btnClose->GetWidgetName() );

    vtkKWText* text = scrolledText->GetWidget();
    text->QuickFormattingOn();
    text->AppendText( "__Mouse Commands in Display View:__\n" );
    text->AppendText( "\n" );

    const char* sGdfButton = 
      QdecUtilities::GetQdecrcResourceString( "GDF_BUTTON" );
    if (sGdfButton && (0 == strcmp(sGdfButton,"BUTTON_1"))) {
      text->AppendText( "**Button 1:** Select vertex (and graph GDF)\n" );
    }
    else {
      text->AppendText( "**Button 1:**                Rotate camera\n" );
    }
    text->AppendText( "**Button 2:**                Pan camera\n" );
    text->AppendText( "**Button 3:**                Zoom camera\n" );
    if (sGdfButton && (0 == strcmp(sGdfButton,"BUTTON_1"))) {
      text->AppendText( "**Ctrl-Button 1:**           Rotate camera\n" );
    }
    else {
      text->AppendText( "**Ctrl-Button 1:** Select vertex (and graph GDF)\n" );
    }
    text->AppendText( "**Shift-Button 1, drag:**    Draw a path\n" );
    text->AppendText( "**Shift-Button 1, no drag:** Clear path\n" );
    text->AppendText( "\n" );
    text->AppendText( "For more help, please visit:\n" );
    text->AppendText("https://surfer.nmr.mgh.harvard.edu/fswiki/Qdec\n");
  }

  mDlogHelp->Display();
}
