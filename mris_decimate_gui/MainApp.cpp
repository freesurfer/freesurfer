/**
 * @file  MainApp.cpp
 * @brief The application that starts all.
 *
 * Main application class for mris_decimate_gui app.  Based off of code
 * from freeview.
 */
/*
 * Original Author: Dan Ginsburg
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:30 $
 *    $Revision: 1.4 $
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

#ifdef __WXMAC__
#include <ApplicationServices/ApplicationServices.h>
#endif
#include "wx/wxprec.h"
#include "MainApp.h"
#include "MainWindow.h"
#include "MyCmdLineParser.h"
#include <wx/wx.h>
#include <wx/xrc/xmlres.h>
#include <wx/config.h>
#include <wx/tooltip.h>
#include <wx/datetime.h>
#include <wx/fs_mem.h>
#include <stdlib.h>
extern "C"
{
#include "error.h"
}

BEGIN_EVENT_TABLE( MainApp, wxApp )
    // EVT_KEY_UP( MainApp::OnKeyDown )
END_EVENT_TABLE()

void MyErrorExit( int nCode )
{
    throw nCode;
}

wxString Progname = _("mris_decimate_gui");

IMPLEMENT_APP( MainApp )

MainApp::MainApp() :
        m_wndMain( NULL )
{
#ifdef __WXMAC__
    ProcessSerialNumber PSN;
    GetCurrentProcess(&PSN);
    TransformProcessType(&PSN,kProcessTransformToForegroundApplication);
#endif
}

MainApp::~MainApp()
{}

bool MainApp::OnInit()
{
    // set error exit
    ErrorSetExitFunc( &MyErrorExit );

    InitializeGUI();

    CmdLineEntry cmdLineDesc[] =
    {
        CmdLineEntry( CMD_LINE_OPTION, "f", "surface", "<FILE>", "Load a surface file on startup.", 1, 1 ),
        CmdLineEntry( CMD_LINE_OPTION, "ss", "screenshot", "<PNG_FILE>", "Take a screen shot (.png) of the viewport and then quit the program.", 1, 1 ),
        CmdLineEntry( CMD_LINE_OPTION, "d", "decimationLevel", "<level>", "Decimation level between (0, 1.0) to apply to surface.", 1, 1 ),
        CmdLineEntry( CMD_LINE_OPTION, "c", "curvature", "<curvatureName>", "Curvature name, can be: (None, Gaussian, Mean, K1, K2, S, C, SI, BE, or FI)", 1, 1),
        CmdLineEntry( CMD_LINE_OPTION, "fs", "filesave", "<FILE>", "Save decimated surface to a file.", 1, 1 ),
        CmdLineEntry( CMD_LINE_OPTION, "r", "rotate", "<angle>", "Rotate (azimuth) the camera <angle> degrees from its initial angle.  Useful for screenshots.", 1, 1),
        CmdLineEntry( CMD_LINE_NONE )
    };

    char progDesc[] = "Surface decimator for freesurfer.";

    MyCmdLineParser cmd( (const char*)"mris_decimate_gui", (CmdLineEntry*)cmdLineDesc );
    cmd.SetProgramDescription( progDesc );
    if ( !cmd.Parse( argc, (char**)argv ) )
    {
        return false;
    }

    string_array sa;

    if ( cmd.Found( "d", &sa ) )
    {
        wxArrayString script;
        script.Add( _( "decimationLevel" ) );
        script.Add( wxString::FromAscii(sa[0].c_str()) );
        m_wndMain->AddScript( script );
    }

    if ( cmd.Found( "f", &sa ) )
    {
        wxArrayString script;
        script.Add( _( "loadsurface" ) );
        script.Add( wxString::FromAscii(sa[0].c_str()) );
        m_wndMain->AddScript( script );
    }

    if ( cmd.Found( "r", &sa ) )
    {
        wxArrayString script;
        script.Add( _( "rotate" ) );
        script.Add( wxString::FromAscii(sa[0].c_str()) );
        m_wndMain->AddScript( script );
    }

    if ( cmd.Found( "c", &sa ) )
    {
        wxArrayString script;
        script.Add( _( "curvature" ) );
        script.Add( wxString::FromAscii(sa[0].c_str()) );
        m_wndMain->AddScript( script );
    }

    if ( cmd.Found( "fs", &sa ) )
    {
        wxArrayString script;
        script.Add( _( "filesave" ) );
        script.Add( wxString::FromAscii(sa[0].c_str()) );
        m_wndMain->AddScript( script );
    }

    if ( cmd.Found( "ss", &sa ) )
    {
        wxArrayString script;
        script.Add( _( "screencapture" ) );
        script.Add( wxString::FromAscii(sa[0].c_str()) );
        m_wndMain->AddScript( script );
        script.Clear();
        script.Add( _( "quit" ) );
        m_wndMain->AddScript( script );
    }

    m_wndMain->Hide();

    m_wndMain->Show();
    m_wndMain->RunScript();
    return true;
}

void MainApp::InitializeGUI()
{
    SetVendorName( _("Children's Hospital Boston") );
    SetAppName( Progname );

    CreateMainWindow();
}

void MainApp::CreateMainWindow()
{
#ifdef __WXMAC__
    wxApp::s_macAboutMenuItemId = ID_ABOUT;
#endif

    // create the main application window
    m_wndMain = new MainWindow( (wxWindow*)NULL );
    m_wndMain->SetTitle( Progname );

    // and show it (the frames, unlike simple controls, are not shown when
    // created initially)
    SetTopWindow( m_wndMain );
    m_wndMain->Show( true );

}

int MainApp::OnExit()
{
    // clean up: Set() returns the active config object as Get() does, but unlike
    // Get() it doesn't try to create one if there is none (definitely not what
    // we want here!)
    delete wxConfigBase::Set( (wxConfigBase *) NULL );

    return 0;
}
