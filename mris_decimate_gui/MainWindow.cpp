/*
 * Original Author: Dan Ginsburg (@ Children's Hospital Boston)
 * CVS Revision Info:
 *    $Author: ginsburg $
 *    $Date: 2010/08/06 14:23:57 $
 *    $Revision: 1.2 $
 *
 * Copyright (C) 2010,
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

#include "MainWindow.h"
#include "DecimatePanel.h"
#include "RenderPanel.h"

#include <wx/panel.h>
#include <wx/filedlg.h>
#include <wx/splitter.h>
#include <wx/config.h>
#include <wx/filename.h>
#include <wx/aboutdlg.h>

static char vcid[] = "$Id: MainWindow.cpp,v 1.2 2010/08/06 14:23:57 ginsburg Exp $";

///////////////////////////////////////////////////////////////////////////////////
//
//  Constructor/Destructor
//
//

///
/// Constructor
///
MainWindow::MainWindow( wxWindow* parent ) :
        MainWindowBase( parent ),
        m_origSurface(NULL),
        m_currentFilePath("")
{

    SetSize( 50, 50, 1024, 768 );

    wxSplitterWindow *splitter = new wxSplitterWindow(this, -1);


    wxBoxSizer *sizer = new wxBoxSizer(wxHORIZONTAL);
    m_renderPanel = new RenderPanel(splitter);
    sizer->Add(m_renderPanel, 1, wxEXPAND);

    m_decimatePanel = new DecimatePanel(splitter, m_renderPanel);

    splitter->SplitVertically(m_decimatePanel, m_renderPanel);
    splitter->SetSashPosition(350);

    wxConfigBase *config = wxConfigBase::Get();
    if (config )
    {
        m_lastLoadDir = config->Read( _("/MainWindow/LastLoadDir"), _("") );
        m_lastSaveDir = config->Read( _("/MainWindow/LastSaveDir"), _("") );
        m_decimatePanel->SetLastSaveDir(&m_lastSaveDir);
    }

#ifdef __WXMAC__
    wxApp:: s_macAboutMenuItemId = ID_ABOUT;
#endif
}

///////////////////////////////////////////////////////////////////////////////////
//
//  wxWidgets Handlers
//
//

void MainWindow::OnFileOpen( wxCommandEvent& event )
{
    if (m_origSurface != NULL)
    {
        MRISfree(&m_origSurface);
        m_origSurface = NULL;
    }

    wxFileDialog* openDialog = new wxFileDialog(this,
            _("Choose a surface file to open"), m_lastLoadDir, wxEmptyString,
            _("All files (*)|*"),
            wxFD_OPEN, wxDefaultPosition);

    // Creates a "open file" dialog
    if (openDialog->ShowModal() == wxID_OK) // if the user click "Open" instead of "Cancel"
    {
        wxString filePath = openDialog->GetPath();
        LoadSurface(filePath);

        wxFileName fileName(filePath);
        m_lastLoadDir = fileName.GetPath();
    }


    // Clean up after ourselves
    openDialog->Destroy();

}

void MainWindow::OnFileExit( wxCommandEvent& event )
{
    Close();
}

void MainWindow::OnFileSave( wxCommandEvent& event )
{
    if (m_currentFilePath == "")
    {
        OnFileSaveAs(event);
        return;
    }

    if (m_decimatePanel->SaveDecimatedSurface(m_currentFilePath.c_str()) != 0)
    {
        wxMessageBox(wxString::Format("ERROR: Saving file '%s'", m_currentFilePath.c_str()), "ERROR");
    }
    wxFileName fileName(m_currentFilePath);
    m_lastSaveDir = fileName.GetPath();
}

void MainWindow::OnFileSaveAs( wxCommandEvent& event )
{
    wxFileDialog *saveDialog = new wxFileDialog(
        this, _("Save File As..."), m_lastSaveDir, wxEmptyString,
        _("All files(*|*"),
        wxFD_SAVE | wxFD_OVERWRITE_PROMPT, wxDefaultPosition);

    if (saveDialog->ShowModal() == wxID_OK)
    {
        m_currentFilePath = saveDialog->GetPath();
        SetTitle(wxString::Format("mris_decimate_gui [%s]", m_currentFilePath.c_str()));
        wxFileName fileName(m_currentFilePath);
        m_lastSaveDir = fileName.GetPath();

        if (m_decimatePanel->SaveDecimatedSurface(m_currentFilePath.c_str()) != 0)
        {
            wxMessageBox(wxString::Format("ERROR: Saving file '%s'", m_currentFilePath.c_str()),
                         "ERROR");
        }
    }
}

void MainWindow::OnAbout( wxCommandEvent& event)
{

    wxAboutDialogInfo info;
    info.SetName("mris_decimate_gui");
    info.SetVersion(vcid);
    info.SetDescription("This program provides tools to decimate and visualize freesurfer surfaces.");
    info.SetCopyright("2010 Children's Hospital Boston - Daniel Ginsburg <daniel.ginsburg@childrens.harvard.edu>");
    wxAboutBox(info);
}

void MainWindow::OnClose( wxCloseEvent& event )
{
    wxConfigBase *config = wxConfigBase::Get();
    if (config)
    {
        config->Write( _("/MainWindow/LastLoadDir"), m_lastLoadDir );
        config->Write( _("/MainWindow/LastSaveDir"), m_lastSaveDir );
    }
    event.Skip();
}

///////////////////////////////////////////////////////////////////////////////////
//
//  Protected Methods
//
//

void MainWindow::LoadSurface(const wxString& filePath)
{
    m_origSurface = MRISfastRead(filePath.c_str());

    if (m_origSurface == NULL)
    {
        wxMessageBox(wxString::Format("ERROR: Loading file '%s'", filePath.c_str()), "ERROR");
    }
    else
    {
        m_decimatePanel->SetOrigSurface(m_origSurface);
        m_decimatePanel->DoDecimate();
        m_renderPanel->ResetCamera();

        m_saveSurface->Enable();
        m_saveSurfaceAs->Enable();
        m_currentFilePath = "";
        SetTitle(wxString::Format("mris_decimate_gui"));
    }
}

///////////////////////////////////////////////////////////////////////////////////
//
//  Public Methods
//
//

void MainWindow::AddScript( const wxArrayString& script )
{
    m_scripts.push_back( script );
}

void MainWindow::RunScript()
{
    if ( m_scripts.size() == 0 )
        return;

    while (!m_scripts.empty() )
    {
        wxArrayString sa = m_scripts[0];
        m_scripts.erase( m_scripts.begin() );
        if ( sa[0] == _("loadsurface") )
        {
            CommandLoadSurface( sa );
        }
        else if ( sa[0] == _("screencapture") )
        {
            CommandScreenCapture( sa );
        }
        else if ( sa[0] == _("decimationLevel") )
        {
            CommandDecimationLevel( sa );
        }
        else if ( sa[0] == _("curvature") )
        {
            CommandCurvature( sa );
        }
        else if ( sa[0] == _("filesave") )
        {
            CommandFileSave( sa );
        }
        else if ( sa[0] == _("rotate") )
        {
            CommandCameraRotate( sa );
        }
        else if ( sa[0] == _("quit") || sa[0] == _("exit") )
        {
            Close();
        }
    }
}

void MainWindow::CommandLoadSurface( const wxArrayString& cmd )
{
    LoadSurface( cmd[1] );

}

void MainWindow::CommandScreenCapture( const wxArrayString& cmd )
{
    m_renderPanel->ShowHistogram( false );
    m_renderPanel->SaveScreenshot( cmd[1] );
}

void MainWindow::CommandDecimationLevel( const wxArrayString& sa )
{
    double dValue;
    if ( sa[1].ToDouble( &dValue ) )
    {
        if ( dValue > 1.0 )
        {
            dValue = 1.0;
        }
        if ( dValue < 0.0 )
        {
            dValue = 0.0;
        }
        m_decimatePanel->UpdateDecimationLevel( (float) dValue );
    }
    else
    {
        std::cerr << "Decimation level is not valid." << endl;
    }
}

void MainWindow::CommandCurvature( const wxArrayString& sa )
{
    if (!m_decimatePanel->UpdateCurvature( sa[1] ) )
    {
        std::cerr << "Unknown curvature type: " << sa[1] << endl;
    }
}

void MainWindow::CommandFileSave( const wxArrayString& sa )
{
    if (m_decimatePanel->SaveDecimatedSurface(sa[1].c_str()) != 0)
    {
        std::cerr << "Error saving decimated surface to file: " << sa[1] << endl;
    }
}

void MainWindow::CommandCameraRotate( const wxArrayString& sa )
{
    double dValue;
    if ( sa[1].ToDouble( &dValue ) )
    {
        m_renderPanel->AzimuthCamera( dValue );
    }
    else
    {
        std::cerr << "Camera rotate angle is not valid." << endl;
    }    
}


