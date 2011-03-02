/**
 * @file  DialogWriteMovieFrames.h
 * @brief Dialog to save volume.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:01 $
 *    $Revision: 1.3 $
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



#include "DialogWriteMovieFrames.h"
#include <wx/xrc/xmlres.h>
#include <wx/filedlg.h>
#include <wx/filename.h>
#include <wx/file.h>
#include <wx/combobox.h>
#include "MainWindow.h"
#include "MyUtils.h"

BEGIN_EVENT_TABLE( DialogWriteMovieFrames, wxDialog )
  EVT_BUTTON    ( XRCID( "ID_BUTTON_CLOSE" ),     DialogWriteMovieFrames::OnClose )
  EVT_BUTTON    ( XRCID( "ID_BUTTON_WRITE" ),     DialogWriteMovieFrames::OnWrite )
  EVT_BUTTON    ( XRCID( "ID_BUTTON_ABORT" ),     DialogWriteMovieFrames::OnAbort )   
  EVT_BUTTON    ( XRCID( "ID_BUTTON_OPEN" ),      DialogWriteMovieFrames::OnOpen )  
  EVT_SHOW      ( DialogWriteMovieFrames::OnShow ) 
END_EVENT_TABLE()


DialogWriteMovieFrames::DialogWriteMovieFrames( wxWindow* parent )
{
  wxXmlResource::Get()->LoadDialog( this, parent, wxT("ID_DIALOG_WRITE_MOVIE_FRAMES") );
  m_textOutputLocation      = XRCCTRL( *this, "ID_TEXT_OUTPUT_LOCATION",  wxTextCtrl );
  m_choiceOutputExtension   = XRCCTRL( *this, "ID_CHOICE_EXTENSION",      wxChoice );
  m_textAngleStep           = XRCCTRL( *this, "ID_TEXT_ANGLE_STEP",       wxTextCtrl );
  m_buttonAbort             = XRCCTRL( *this, "ID_BUTTON_ABORT",          wxButton );
  m_buttonWrite             = XRCCTRL( *this, "ID_BUTTON_WRITE",          wxButton );
  m_buttonClose             = XRCCTRL( *this, "ID_BUTTON_CLOSE",          wxButton );
}

DialogWriteMovieFrames::~DialogWriteMovieFrames()
{}

wxString DialogWriteMovieFrames::GetOutputLocation()
{
  return m_textOutputLocation->GetValue().Trim( true ).Trim( false );
}

wxString DialogWriteMovieFrames::GetOutputExtension()
{
  return m_choiceOutputExtension->GetStringSelection();
}

double DialogWriteMovieFrames::GetAngleStep()
{
  double dvalue;
  if ( m_textAngleStep->GetValue().ToDouble( &dvalue ) )
    return dvalue;
  else
    return 0;
}

void DialogWriteMovieFrames::OnClose( wxCommandEvent& event )
{
  Close();
}

void DialogWriteMovieFrames::OnWrite( wxCommandEvent& event )
{
  if ( GetOutputLocation().IsEmpty() )
  {
    wxMessageDialog dlg( this, 
                         _("Output location can not be empty."), 
                         _("Error"), wxOK | wxICON_ERROR );
    dlg.ShowModal();
    return;
  }
  if ( MainWindow::GetMainWindowPointer()->GetMainViewId() == 3 && GetAngleStep() == 0 )
  {
    wxMessageDialog dlg( this, 
                         _("Angle step can not be 0."), 
                         _("Error"), wxOK | wxICON_ERROR );
    dlg.ShowModal();
    return;
  }
  
  MainWindow::GetMainWindowPointer()->StartWriteMovieFrames();
  UpdateUI();
}

void DialogWriteMovieFrames::OnAbort( wxCommandEvent& event )
{
  MainWindow::GetMainWindowPointer()->StopWriteMovieFrames();
  UpdateUI();
}

void DialogWriteMovieFrames::OnOpen( wxCommandEvent& event )
{
  wxDirDialog dlg( this, 
                   _("Choose a directory"), 
                   GetOutputLocation() );
  if ( dlg.ShowModal() == wxID_OK )
  {
    m_textOutputLocation->SetValue( dlg.GetPath() );
  }
}

void DialogWriteMovieFrames::OnShow( wxShowEvent& event )
{
  UpdateUI(); 
}

void DialogWriteMovieFrames::UpdateUI()
{
  bool bWriting = MainWindow::GetMainWindowPointer()->IsWritingMovieFrames();
  m_buttonAbort->Enable( bWriting );
  m_buttonWrite->Enable( !bWriting );
  m_buttonClose->Enable( !bWriting );
  m_textAngleStep->Enable( MainWindow::GetMainWindowPointer()->GetMainViewId() == 3 );
}
