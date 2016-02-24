/*
 * Original Author: Dan Ginsburg (@ Children's Hospital Boston)
 * CVS Revision Info:
 *    $Author: zkaufman $
 *    $Date: 2016/02/24 16:28:03 $
 *    $Revision: 1.5 $
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

#include "MainWindow.h"
#include "DecimatePanel.h"
#include "RenderPanel.h"

#include <wx/panel.h>
#include <wx/filedlg.h>
#include <wx/splitter.h>
#include <wx/config.h>
#include <wx/filename.h>
#include <wx/aboutdlg.h>

///////////////////////////////////////////////////////////////////////////////////
//
//  VTK Observer classes
//
//


class CameraObserver : public vtkCommand
{
public:
    static CameraObserver *New()
    {
        return new CameraObserver;
    }
 
    void SetMainWindow( MainWindow *mainWindow )
    {
        m_mainWindow = mainWindow;
    }

    virtual void Execute(vtkObject *vtkNotUsed(caller),
                         unsigned long event,
                         void *vtkNotUsed(calldata))
    {
        if (m_mainWindow != 0)
        {
            m_mainWindow->HandleCameraUpdate();
        }
    }
 
protected:
    CameraObserver() :
        m_mainWindow(0)
    {
    }

    ~CameraObserver()
    {
        m_mainWindow = 0;
    }

    MainWindow *m_mainWindow;
};

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
        m_currentFilePath(_(""))
{
    SetSize( 50, 50, 1024, 840 );

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

    // Add camera observer
    vtkSmartPointer<CameraObserver> cameraObserver = 
        vtkSmartPointer<CameraObserver>::New();
    cameraObserver->SetMainWindow(this);
    m_renderPanel->AddCameraObserver(cameraObserver);
    

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
    if (m_currentFilePath == _(""))
    {
        OnFileSaveAs(event);
        return;
    }

    if (m_decimatePanel->SaveDecimatedSurface(m_currentFilePath.mb_str(wxConvUTF8)) != 0)
    {
        wxMessageBox(wxString::Format(_("ERROR: Saving file '%s'"), m_currentFilePath.c_str()), _("ERROR"));
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
        SetTitle(wxString::Format(_("mris_decimate_gui [%s]"), m_currentFilePath.c_str()));
        wxFileName fileName(m_currentFilePath);
        m_lastSaveDir = fileName.GetPath();

        if (m_decimatePanel->SaveDecimatedSurface(m_currentFilePath.mb_str(wxConvUTF8)) != 0)
        {
            wxMessageBox(wxString::Format(_("ERROR: Saving file '%s'"), m_currentFilePath.c_str()),
                         _("ERROR"));
        }
    }
}

void MainWindow::OnAbout( wxCommandEvent& event)
{

    wxAboutDialogInfo info;
    info.SetName(_("mris_decimate_gui"));
	wxString version =
			wxString::FromAscii( __DATE__) +
    		_(" ") +
		    wxString::FromAscii(__TIME__);

    info.SetVersion( _("1.0 (internal) \r\nbuild ") + version);
    info.SetDescription(_("This program provides tools to decimate and visualize freesurfer surfaces."));
    info.SetCopyright(_("2010 Children's Hospital Boston\n"
	    			  "Daniel Ginsburg <daniel.ginsburg@childrens.harvard.edu>\n"
	    			  "Rudolph Pienaar <rudolph.pienaar@childrens.harvard.edu>"));
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
    m_origSurface = MRISfastRead(filePath.mb_str(wxConvUTF8));

    if (m_origSurface == NULL)
    {
        wxMessageBox(wxString::Format(_("ERROR: Loading file '%s'"), filePath.c_str()), _("ERROR"));
    }
    else
    {
        m_decimatePanel->SetOrigSurface(m_origSurface);
        m_decimatePanel->DoDecimate();
        m_renderPanel->ResetCamera();
        m_renderPanel->RecordCameraCoordinates();

        m_saveSurface->Enable();
        m_saveSurfaceAs->Enable();
        m_currentFilePath = _("");
        SetTitle(wxString::Format(_("mris_decimate_gui")));
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

/// Handle updates to the camera
void MainWindow::HandleCameraUpdate()
{
    m_decimatePanel->UpdateCameraInfo();
}

void MainWindow::CommandLoadSurface( const wxArrayString& cmd )
{
    LoadSurface( cmd[1] );

}

void MainWindow::CommandScreenCapture( const wxArrayString& cmd )
{
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
    if (m_decimatePanel->SaveDecimatedSurface(sa[1].mb_str(wxConvUTF8)) != 0)
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


