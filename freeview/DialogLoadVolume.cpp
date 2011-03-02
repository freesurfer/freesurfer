/**
 * @file  DialogLoadVolume.h
 * @brief Dialog to load volume data.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:01 $
 *    $Revision: 1.21 $
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



#include "DialogLoadVolume.h"
#include <wx/xrc/xmlres.h>
#include <wx/filedlg.h>
#include <wx/filename.h>
#include "MyUtils.h"
#include "MainWindow.h"
#include "LayerMRI.h"
#include "LayerPropertiesMRI.h"
#include "LUTDataHolder.h"

BEGIN_EVENT_TABLE( DialogLoadVolume, wxDialog )
  EVT_BUTTON    ( wxID_OK,                        DialogLoadVolume::OnOK )
  EVT_BUTTON    ( XRCID( "ID_BUTTON_FILE" ),      DialogLoadVolume::OnButtonOpen )
  EVT_COMBOBOX  ( XRCID( "ID_COMBO_FILENAME" ),   DialogLoadVolume::OnFileSelectionChanged )
  EVT_BUTTON    ( XRCID( "ID_BUTTON_REG_FILE" ),  DialogLoadVolume::OnButtonRegFile )
  EVT_CHECKBOX  ( XRCID( "ID_CHECK_APPLY_REG" ),  DialogLoadVolume::OnCheckApplyReg )
  EVT_CHOICE    ( XRCID( "ID_CHOICE_COLORMAP" ),  DialogLoadVolume::OnChoiceColorMap )
  EVT_CHOICE    ( XRCID( "ID_CHOICE_LUT" ),       DialogLoadVolume::OnChoiceLUT )
END_EVENT_TABLE()


DialogLoadVolume::DialogLoadVolume( wxWindow* parent, bool bEnableResample )
{
  wxXmlResource::Get()->LoadDialog( this, parent, 
				    wxT("ID_DIALOG_LOAD_VOLUME") );
  m_checkResample   = XRCCTRL( *this, "ID_CHECK_RESAMPLE",  wxCheckBox );
  m_btnOpen         = XRCCTRL( *this, "ID_BUTTON_FILE",     wxButton );
  m_comboFileName   = XRCCTRL( *this, "ID_COMBO_FILENAME",  wxComboBox );
  m_checkApplyReg   = XRCCTRL( *this, "ID_CHECK_APPLY_REG", wxCheckBox );
  m_textRegFile     = XRCCTRL( *this, "ID_TEXT_REG_FILE",   wxTextCtrl );
  m_btnRegFile      = XRCCTRL( *this, "ID_BUTTON_REG_FILE", wxButton );
  m_radioNearest    = XRCCTRL( *this, "ID_RADIO_NEAREST",   wxRadioButton );
  m_radioTrilinear  = XRCCTRL( *this, "ID_RADIO_TRILINEAR", wxRadioButton );
  m_choiceColorMap  = XRCCTRL( *this, "ID_CHOICE_COLORMAP", wxChoice );
  m_choiceLUT       = XRCCTRL( *this, "ID_CHOICE_LUT",      wxChoice );
  m_staticLUT       = XRCCTRL( *this, "ID_STATIC_LUT",      wxStaticText );
  m_checkResample->Show( bEnableResample );
  m_comboFileName->SetFocus();
  UpdateLUT();
}

DialogLoadVolume::~DialogLoadVolume()
{}

void DialogLoadVolume::UpdateLUT()
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

wxArrayString DialogLoadVolume::GetVolumeFileNames()
{
  return MyUtils::SplitString( m_comboFileName->GetValue(), _(";") );
}

wxString DialogLoadVolume::GetRegFileName()
{
  if ( m_checkApplyReg->IsChecked() )
    return m_textRegFile->GetValue().Trim( true ).Trim( false );
  else
    return _("");
}

void DialogLoadVolume::OnOK( wxCommandEvent& event )
{
  if ( GetVolumeFileNames().IsEmpty())
  {
    wxMessageDialog dlg( this, 
			 _("Volume file names cannot be empty."), 
			 _("Error"), 
			 wxOK | wxICON_ERROR );
    dlg.ShowModal();
    return;
  }
  else if ( m_checkApplyReg->IsChecked() && 
	    m_textRegFile->GetValue().Trim( true ).Trim( false ).IsEmpty() )
  {
    wxMessageDialog dlg( this, 
			 _("Registration file name cannot be empty."), 
			 _("Error"), 
			 wxOK | wxICON_ERROR );
    dlg.ShowModal();
    return;
  }

  event.Skip();
}

bool DialogLoadVolume::IsToResample()
{
  return m_checkResample->IsChecked();
}

int DialogLoadVolume::GetSampleMethod()
{
  if ( m_radioNearest->GetValue() )
    return 0;         
  else
    return 1;
} 

void DialogLoadVolume::SetRecentFiles( const wxArrayString& list )
{
  m_comboFileName->Clear();
  for ( int i = 0; i < (int)list.GetCount(); i++ )
  {
    m_comboFileName->Append( list[i] );
  }
  if ( list.GetCount() > 0 )
  {
    m_comboFileName->SetSelection( 0 );
    m_comboFileName->SetInsertionPointEnd();
  }
}

void DialogLoadVolume::OnButtonOpen( wxCommandEvent& event )
{
  wxFileDialog dlg
    ( this, 
      _("Open volume file"), 
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
    m_comboFileName->SetValue( text );
    m_comboFileName->SetInsertionPointEnd();
    m_strLastDir = wxFileName( dlg.GetPath() ).GetPath();
  }
}

void DialogLoadVolume::OnButtonRegFile( wxCommandEvent& event )
{
  wxArrayString fns = GetVolumeFileNames();
  if ( !fns.IsEmpty() )
    m_strLastDir = wxFileName( fns[0] ).GetPath();
  wxFileDialog dlg
    ( this, 
      _("Open registration file"), 
      m_strLastDir, _(""),
      _("Registration files (*.dat;*.xfm;*.lta;*.mat)|*.dat;*.xfm;*.lta;*.mat|All files (*.*)|*.*"),
      wxFD_OPEN );
  if ( dlg.ShowModal() == wxID_OK )
  {
    m_textRegFile->SetValue( dlg.GetPath() );
    m_textRegFile->SetInsertionPointEnd();
    // m_strLastDir = wxFileName( dlg.GetPath() ).GetPath();
  }
}

void DialogLoadVolume::OnFileSelectionChanged( wxCommandEvent& event )
{
  m_comboFileName->SetInsertionPointEnd();
}

void DialogLoadVolume::OnCheckApplyReg( wxCommandEvent& event )
{
  m_textRegFile->Enable( event.IsChecked() );
  m_btnRegFile->Enable( event.IsChecked() );
}

void DialogLoadVolume::OnChoiceColorMap( wxCommandEvent& event )
{
  m_staticLUT->Enable( event.GetSelection() == LayerPropertiesMRI::LUT );
  m_choiceLUT->Enable( event.GetSelection() == LayerPropertiesMRI::LUT );
}

void DialogLoadVolume::OnChoiceLUT( wxCommandEvent& event )
{
  LUTDataHolder* luts = MainWindow::GetMainWindowPointer()->GetLUTData();
  if ( event.GetSelection() >= luts->GetCount() )
  {
    wxFileDialog dlg( this, _("Load lookup table file"), m_strLastDir, _(""),
                      _("LUT files (*.*)|*.*"),
                      wxFD_OPEN );
    if ( dlg.ShowModal() == wxID_OK && luts->LoadColorTable( dlg.GetPath().c_str() ) )
    {
      UpdateLUT();
      m_choiceLUT->SetSelection( luts->GetCount() - 1 );
    }
  }
}

wxString DialogLoadVolume::GetColorMap()
{
  const char* names[] = { "grayscale", "lut", "heat", "jet", "gecolor", "nih" };
  return names[m_choiceColorMap->GetSelection()];
}

wxString DialogLoadVolume::GetLUT()
{
  return m_choiceLUT->GetString( m_choiceLUT->GetSelection() );
}
