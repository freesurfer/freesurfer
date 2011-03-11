/**
 * @file  DialogLoadPVolumes.h
 * @brief Dialog to load DTI data.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: krish $
 *    $Date: 2011/03/11 23:27:36 $
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



#include "DialogLoadPVolumes.h"
#include <wx/xrc/xmlres.h>
#include <wx/filedlg.h>
#include <wx/filename.h>
#include "MainWindow.h"
#include "LUTDataHolder.h"
#include "MyUtils.h"

BEGIN_EVENT_TABLE( DialogLoadPVolumes, wxDialog )
  EVT_BUTTON    ( wxID_OK,                          DialogLoadPVolumes::OnOK )
  EVT_BUTTON    ( XRCID( "ID_BUTTON_OPEN" ),        DialogLoadPVolumes::OnButtonOpen )
  EVT_CHOICE    ( XRCID( "ID_CHOICE_LUT" ),         DialogLoadPVolumes::OnChoiceLUT )
END_EVENT_TABLE()


DialogLoadPVolumes::DialogLoadPVolumes( wxWindow* parent )
{
  wxXmlResource::Get()->LoadDialog( this, parent, wxT("ID_DIALOG_LOAD_PVOLUMES") );
  m_textFileNames   = XRCCTRL( *this, "ID_TEXT_VOLUME_FILES",     wxTextCtrl );
  m_textPrefix      = XRCCTRL( *this, "ID_TEXT_FILENAME_PREFIX",  wxTextCtrl );
  m_choiceLUT       = XRCCTRL( *this, "ID_CHOICE_LUT",  wxChoice );
  m_btnOpen         = XRCCTRL( *this, "ID_BUTTON_OPEN", wxButton );
  m_textFileNames->SetFocus();
  
  UpdateLUT();
}

DialogLoadPVolumes::~DialogLoadPVolumes()
{}


void DialogLoadPVolumes::UpdateLUT()
{
  LUTDataHolder* luts = MainWindow::GetMainWindowPointer()->GetLUTData();
  m_choiceLUT->Clear();
  for ( int i = 0; i < luts->GetCount(); i++ )
  {
    m_choiceLUT->Append( wxString::FromAscii( luts->GetName( i ) ) );
  }
  m_choiceLUT->Append( _("Load lookup table...") );
  m_choiceLUT->SetSelection( 0 );
}

void DialogLoadPVolumes::OnChoiceLUT( wxCommandEvent& event )
{
  LUTDataHolder* luts = MainWindow::GetMainWindowPointer()->GetLUTData();
  if ( event.GetSelection() >= luts->GetCount() )
  {
    wxFileDialog dlg( this, _("Load lookup table file"), _(""), _(""),
                      _("LUT files (*.*)|*.*"),
                      wxFD_OPEN );
    if ( dlg.ShowModal() == wxID_OK && luts->LoadColorTable( dlg.GetPath().c_str() ) )
    {
      UpdateLUT();
      m_choiceLUT->SetSelection( luts->GetCount() - 1 );
    }
  }
}

wxArrayString DialogLoadPVolumes::GetVolumeFileNames()
{
  return MyUtils::SplitString( m_textFileNames->GetValue(), _(";") );
}

wxString DialogLoadPVolumes::GetLUT()
{
  return m_choiceLUT->GetString( m_choiceLUT->GetSelection() );
}

wxString DialogLoadPVolumes::GetFileNamePrefix()
{
  return m_textPrefix->GetValue().Trim( true ).Trim( false );
}

void DialogLoadPVolumes::OnOK( wxCommandEvent& event )
{
  if ( GetVolumeFileNames().Count() == 0 )
  {
    wxMessageDialog dlg
      ( this, 
	_("Please select volumes to load."), 
	_("Error"), wxOK | wxICON_ERROR );
    dlg.ShowModal();
    return;
  }
  
  event.Skip();
}


void DialogLoadPVolumes::OnButtonOpen( wxCommandEvent& event )
{
  wxFileDialog dlg
      ( this, 
        _("Open p-label volume files"), 
        m_strLastDir, _(""),
        _("Volume files (*.nii;*.nii.gz;*.img;*.mgz)|*.nii;*.nii.gz;*.img;*.mgz|All files (*.*)|*.*"),
        wxFD_OPEN | wxFD_MULTIPLE );
  if ( dlg.ShowModal() == wxID_OK )
  {
    wxArrayString fns;
    dlg.GetPaths( fns );
    wxString text;
    for ( size_t i = 0; i < fns.GetCount(); i++ )
    {
      text += fns[i];
      if ( i != fns.GetCount()-1 )
        text += _(";");
    }
    m_textFileNames->SetValue( text );
    m_textFileNames->SetInsertionPointEnd();
    m_strLastDir = wxFileName( dlg.GetPath() ).GetPath();
  }
}

