/**
 * @file  DialogSavePointSetAs.h
 * @brief Dialog to save point set.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: krish $
 *    $Date: 2011/03/11 23:27:37 $
 *    $Revision: 1.1 $
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



#include "DialogSavePointSetAs.h"
#include <wx/xrc/xmlres.h>
#include <wx/filedlg.h>
#include <wx/filename.h>
#include <wx/file.h>
#include "LayerPropertiesWayPoints.h"

BEGIN_EVENT_TABLE( DialogSavePointSetAs, wxDialog )
  EVT_BUTTON    ( wxID_OK, DialogSavePointSetAs::OnOK )
  EVT_BUTTON    ( XRCID( "ID_BUTTON_FILENAME" ),    DialogSavePointSetAs::OnButtonOpen )
END_EVENT_TABLE()


DialogSavePointSetAs::DialogSavePointSetAs( wxWindow* parent )
{
  wxXmlResource::Get()->LoadDialog( this, parent, wxT("ID_DIALOG_SAVE_POINTSET_AS") );
  m_textFilename    = XRCCTRL( *this, "ID_TEXT_FILENAME", wxTextCtrl );
  m_btnOpen         = XRCCTRL( *this, "ID_BUTTON_FILENAME", wxButton );
  m_radioWayPoints  = XRCCTRL( *this, "ID_RADIO_WAY_POINTS", wxRadioButton );
  m_radioControlPoints  = XRCCTRL( *this, "ID_RADIO_CONTROL_POINTS", wxRadioButton );
  m_textFilename->SetFocus();
}

DialogSavePointSetAs::~DialogSavePointSetAs()
{}

wxString DialogSavePointSetAs::GetFileName()
{
  return m_textFilename->GetValue().Trim( true ).Trim( false );
}

void DialogSavePointSetAs::SetFileName( const wxString& filename )
{
  m_textFilename->ChangeValue( filename );
}

void DialogSavePointSetAs::OnOK( wxCommandEvent& event )
{
  if ( GetFileName().IsEmpty() )
  {
    wxMessageDialog dlg( this, 
                         _("File name can not be empty."), 
                         _("Error"), wxOK | wxICON_ERROR );
    dlg.ShowModal();
    return;
  }
  else 
  {
    if ( wxFile::Exists( GetFileName() ) )
    {
      
      wxMessageDialog dlg( this, 
                           _("File exists. Do you want to overwrite it?"), 
                           _("Warning"), 
                           wxYES_NO | wxNO_DEFAULT );
      if ( dlg.ShowModal() == wxID_NO )
        return;
    }
  }

  event.Skip();
}

void DialogSavePointSetAs::OnButtonOpen( wxCommandEvent& event )
{
  wxFileDialog dlg( this, _("Select file to save"), 
                    wxFileName( GetFileName() ).GetPath(), 
                    _(""),
                    _("All files (*.*)|*.*"),
                    wxFD_OPEN );
  if ( dlg.ShowModal() == wxID_OK )
  {
    m_textFilename->ChangeValue( dlg.GetPath() );
    m_textFilename->ShowPosition( m_textFilename->GetLastPosition() );
  }
}

void DialogSavePointSetAs::SetPointSetType( int nType )
{
  m_radioWayPoints->SetValue( nType == LayerPropertiesWayPoints::WayPoints );
  m_radioControlPoints->SetValue( nType == LayerPropertiesWayPoints::ControlPoints );
}

int DialogSavePointSetAs::GetPointSetType()
{
  if ( m_radioControlPoints->GetValue() )
    return LayerPropertiesWayPoints::ControlPoints;
  else
    return LayerPropertiesWayPoints::WayPoints;
}
