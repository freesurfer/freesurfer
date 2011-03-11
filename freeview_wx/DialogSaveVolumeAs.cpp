/**
 * @file  DialogSaveVolumeAs.h
 * @brief Dialog to save volume.
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



#include "DialogSaveVolumeAs.h"
#include <wx/xrc/xmlres.h>
#include <wx/filedlg.h>
#include <wx/filename.h>
#include <wx/file.h>

BEGIN_EVENT_TABLE( DialogSaveVolumeAs, wxDialog )
  EVT_BUTTON    ( wxID_OK, DialogSaveVolumeAs::OnOK )
  EVT_BUTTON    ( XRCID( "ID_BUTTON_FILENAME" ),    DialogSaveVolumeAs::OnButtonOpen )
END_EVENT_TABLE()


DialogSaveVolumeAs::DialogSaveVolumeAs( wxWindow* parent )
{
  wxXmlResource::Get()->LoadDialog( this, parent, wxT("ID_DIALOG_SAVE_VOLUME_AS") );
  m_textFilename  = XRCCTRL( *this, "ID_TEXT_FILENAME", wxTextCtrl );
  m_btnOpen       = XRCCTRL( *this, "ID_BUTTON_FILENAME", wxButton );
  m_checkReorient = XRCCTRL( *this, "ID_CHECK_REORIENT", wxCheckBox );
  m_textFilename->SetFocus();
}

DialogSaveVolumeAs::~DialogSaveVolumeAs()
{}

wxString DialogSaveVolumeAs::GetFileName()
{
  return m_textFilename->GetValue().Trim( true ).Trim( false );
}

void DialogSaveVolumeAs::SetFileName( const wxString& filename )
{
  m_textFilename->ChangeValue( filename );
}

void DialogSaveVolumeAs::OnOK( wxCommandEvent& event )
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

void DialogSaveVolumeAs::OnButtonOpen( wxCommandEvent& event )
{
  wxFileDialog dlg( this, _("Select file to save"), 
                    wxFileName( GetFileName() ).GetPath(), 
                    _(""),
                    _("Volume files (*.mgz;*.mgh;*.nii;*.nii.gz;*.img)|*.mgz;*.mgh;*.nii;*.nii.gz;*.img|All files (*.*)|*.*"),
                    wxFD_OPEN );
  if ( dlg.ShowModal() == wxID_OK )
  {
    m_textFilename->ChangeValue( dlg.GetPath() );
    m_textFilename->ShowPosition( m_textFilename->GetLastPosition() );
  }
}


bool DialogSaveVolumeAs::GetReorient()
{
  return m_checkReorient->IsChecked();
}
