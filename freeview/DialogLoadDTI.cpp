/**
 * @file  DialogLoadDTI.h
 * @brief Preferences dialog.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2008/10/08 19:14:35 $
 *    $Revision: 1.3 $
 *
 * Copyright (C) 2002-2008,
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
 


#include "DialogLoadDTI.h"
#include <wx/xrc/xmlres.h>
#include <wx/filedlg.h>
#include <wx/filename.h>

BEGIN_EVENT_TABLE( DialogLoadDTI, wxDialog )
	EVT_BUTTON			( wxID_OK,			 						DialogLoadDTI::OnOK )
	EVT_BUTTON			( XRCID( wxT( "ID_BUTTON_VECTOR_FILE" ) ),	DialogLoadDTI::OnButtonVector )
	EVT_BUTTON			( XRCID( wxT( "ID_BUTTON_FA_FILE" ) ),		DialogLoadDTI::OnButtonFA )
	EVT_COMBOBOX		( XRCID( wxT( "ID_COMBO_FA_FILE" ) ), 		DialogLoadDTI::OnComboFASelectionChanged )
END_EVENT_TABLE()


DialogLoadDTI::DialogLoadDTI( wxWindow* parent ) 
{
	wxXmlResource::Get()->LoadDialog( this, parent, wxT("ID_DIALOG_LOAD_DTI") );	
	m_textVector = XRCCTRL( *this, "ID_TEXT_VECTOR_FILE", wxTextCtrl );
	m_comboFA = XRCCTRL( *this, "ID_COMBO_FA_FILE", wxComboBox );
	m_btnVector = XRCCTRL( *this, "ID_BUTTON_VECTOR_FILE", wxButton );
	m_btnFA = XRCCTRL( *this, "ID_BUTTON_FA_FILE", wxButton );
	m_checkNoResample = XRCCTRL( *this, "ID_CHECK_NO_RESAMPLE", wxCheckBox );
	m_textVector->SetFocus();
}

DialogLoadDTI::~DialogLoadDTI()
{
}
	
wxString DialogLoadDTI::GetVectorFileName()
{
	return m_textVector->GetValue().Trim( true ).Trim( false );
}
	
wxString DialogLoadDTI::GetFAFileName()
{
	return m_comboFA->GetValue().Trim( true ).Trim( false );
}

void DialogLoadDTI::OnOK( wxCommandEvent& event )
{
	if ( GetVectorFileName().IsEmpty() )
	{
		wxMessageDialog dlg( this, "Vector file name can not be empty.", "Error", wxOK | wxICON_ERROR );
		dlg.ShowModal();
		return;
	}
	else if ( GetFAFileName().IsEmpty() )
	{
		wxMessageDialog dlg( this, "FA file name can not be empty.", "Error", wxOK | wxICON_ERROR );
		dlg.ShowModal();
		return;
	}
		
	event.Skip();
}

void DialogLoadDTI::OnButtonVector( wxCommandEvent& event )
{	
	wxFileDialog dlg( this, _("Select vector file"), m_strLastDir, _(""), 
				_T("Volume files (*.nii;*.nii.gz;*.img;*.mgz)|*.nii;*.nii.gz;*.img;*.mgz|All files (*.*)|*.*"), 
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
					  _T("Volume files (*.nii;*.nii.gz;*.img;*.mgz)|*.nii;*.nii.gz;*.img;*.mgz|All files (*.*)|*.*"), 
					  wxFD_OPEN );
	if ( dlg.ShowModal() == wxID_OK )
	{
		m_comboFA->SetValue( dlg.GetPath() );
		m_comboFA->SetInsertionPointEnd();
		m_strLastDir = wxFileName( dlg.GetPath() ).GetPath();
	}
}

bool DialogLoadDTI::IsToResample()
{
	return !m_checkNoResample->IsChecked();
}


void DialogLoadDTI::Initialize( bool bResample, bool bEnableCheckBox )
{
	m_checkNoResample->SetValue( !bResample );
	m_checkNoResample->Enable( bEnableCheckBox );
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
