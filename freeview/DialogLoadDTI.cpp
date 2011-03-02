/**
 * @file  DialogLoadDTI.h
 * @brief Dialog to load DTI data.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:01 $
 *    $Revision: 1.11 $
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



#include "DialogLoadDTI.h"
#include <wx/xrc/xmlres.h>
#include <wx/filedlg.h>
#include <wx/filename.h>

BEGIN_EVENT_TABLE( DialogLoadDTI, wxDialog )
  EVT_BUTTON    ( wxID_OK, DialogLoadDTI::OnOK )
  EVT_BUTTON    ( XRCID( "ID_BUTTON_VECTOR_FILE" ), DialogLoadDTI::OnButtonVector )
  EVT_BUTTON    ( XRCID( "ID_BUTTON_FA_FILE" ),     DialogLoadDTI::OnButtonFA )
  EVT_BUTTON    ( XRCID( "ID_BUTTON_REG_FILE" ),    DialogLoadDTI::OnButtonReg )
  EVT_COMBOBOX  ( XRCID( "ID_COMBO_FA_FILE" ),      DialogLoadDTI::OnComboFASelectionChanged )
  EVT_CHECKBOX  ( XRCID( "ID_CHECK_REG" ),          DialogLoadDTI::OnCheckApplyReg )
END_EVENT_TABLE()


DialogLoadDTI::DialogLoadDTI( wxWindow* parent )
{
  wxXmlResource::Get()->LoadDialog( this, parent, wxT("ID_DIALOG_LOAD_DTI") );
  m_textVector  = XRCCTRL( *this, "ID_TEXT_VECTOR_FILE", wxTextCtrl );
  m_textReg     = XRCCTRL( *this, "ID_TEXT_REG_FILE", wxTextCtrl );
  m_comboFA     = XRCCTRL( *this, "ID_COMBO_FA_FILE", wxComboBox );
  m_btnVector   = XRCCTRL( *this, "ID_BUTTON_VECTOR_FILE", wxButton );
  m_btnFA       = XRCCTRL( *this, "ID_BUTTON_FA_FILE", wxButton );
  m_btnReg      = XRCCTRL( *this, "ID_BUTTON_REG_FILE", wxButton );
  m_checkResample = XRCCTRL( *this, "ID_CHECK_RESAMPLE", wxCheckBox );
  m_checkReg    = XRCCTRL( *this, "ID_CHECK_REG", wxCheckBox );
  m_textVector->SetFocus();
  m_btnReg->Enable( false );
  m_textReg->Enable( false );
}

DialogLoadDTI::~DialogLoadDTI()
{}

wxString DialogLoadDTI::GetVectorFileName()
{
  return m_textVector->GetValue().Trim( true ).Trim( false );
}

wxString DialogLoadDTI::GetFAFileName()
{
  return m_comboFA->GetValue().Trim( true ).Trim( false );
}

wxString DialogLoadDTI::GetRegFileName()
{
  if ( m_checkReg->IsChecked() )
    return m_textReg->GetValue().Trim( true ).Trim( false );
  else
    return _("");
}

void DialogLoadDTI::OnOK( wxCommandEvent& event )
{
  if ( GetVectorFileName().IsEmpty() )
  {
    wxMessageDialog dlg
      ( this, 
	_("Vector file name can not be empty."), 
	_("Error"), wxOK | wxICON_ERROR );
    dlg.ShowModal();
    return;
  }
  else if ( GetFAFileName().IsEmpty() )
  {
    wxMessageDialog dlg
      ( this, 
	_("FA file name can not be empty."), 
	_("Error"), 
	wxOK | wxICON_ERROR );
    dlg.ShowModal();
    return;
  }
  else if ( m_checkReg->IsChecked() && GetRegFileName().IsEmpty() )
  {
    wxMessageDialog dlg
      ( this, 
	_("Registration file name can not be empty."), 
	_("Error"), 
	wxOK | wxICON_ERROR );
    dlg.ShowModal();
    return;
  }

  event.Skip();
}

void DialogLoadDTI::OnButtonVector( wxCommandEvent& event )
{
  wxFileDialog dlg
    ( this, 
      _("Select vector file"), 
      m_strLastDir, 
      _(""),
      _("Volume files (*.nii;*.nii.gz;*.img;*.mgz)|*.nii;*.nii.gz;*.img;*.mgz|All files (*.*)|*.*"),
      wxFD_OPEN );
  if ( dlg.ShowModal() == wxID_OK )
  {
    m_textVector->ChangeValue( dlg.GetPath() );
    m_textVector->SetInsertionPointEnd();
    m_textVector->ShowPosition( m_textVector->GetLastPosition() );
    m_strLastDir = wxFileName( dlg.GetPath() ).GetPath();
  }
}

void DialogLoadDTI::OnButtonFA( wxCommandEvent& event )
{
  wxFileDialog dlg( this, _("Select FA file"), m_strLastDir, _(""),
                    _("Volume files (*.nii;*.nii.gz;*.img;*.mgz)|*.nii;*.nii.gz;*.img;*.mgz|All files (*.*)|*.*"),
                    wxFD_OPEN );
  if ( dlg.ShowModal() == wxID_OK )
  {
    m_comboFA->SetValue( dlg.GetPath() );
    m_comboFA->SetInsertionPointEnd();
    m_strLastDir = wxFileName( dlg.GetPath() ).GetPath();
  }
}

void DialogLoadDTI::OnButtonReg( wxCommandEvent& event )
{
  wxFileDialog dlg( this, _("Select Registration file"), m_strLastDir, _(""),
                    _("Registration files (*.dat;*.xfm;*.lta;*.mat)|*.dat;*.xfm;*.lta;*.mat|All files (*.*)|*.*"),
                    wxFD_OPEN );
  if ( dlg.ShowModal() == wxID_OK )
  {
    m_textReg->ChangeValue( dlg.GetPath() );
    m_textReg->SetInsertionPointEnd();
    m_textReg->ShowPosition( m_textReg->GetLastPosition() );
    // m_strLastDir = wxFileName( dlg.GetPath() ).GetPath();
  }
}


bool DialogLoadDTI::IsToResample()
{
  return m_checkResample->IsChecked();
}


void DialogLoadDTI::Initialize( bool bResample, bool bEnableCheckBox )
{
  m_checkResample->SetValue( bResample );
  m_checkResample->Enable( bEnableCheckBox );
}

void DialogLoadDTI::SetRecentFiles( const wxArrayString& list )
{
  m_comboFA->Clear();
  for ( int i = 0; i < (int)list.GetCount(); i++ )
  {
    m_comboFA->Append( list[i] );
  }
  if ( list.GetCount() > 0 )
  {
    m_comboFA->SetSelection( 0 );
    m_comboFA->SetInsertionPointEnd();
  }
}

void DialogLoadDTI::OnComboFASelectionChanged( wxCommandEvent& event )
{
  wxString strg = wxFileName( GetFAFileName() ).GetPath();
  if ( !strg.IsEmpty() )
    m_strLastDir = strg;
}


void DialogLoadDTI::OnCheckApplyReg( wxCommandEvent& event )
{
  m_textReg->Enable( event.IsChecked() );
  m_btnReg->Enable( event.IsChecked() );
}
