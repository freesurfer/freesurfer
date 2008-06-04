/**
 * @file  DialogLoadVolume.h
 * @brief Preferences dialog.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2008/06/04 20:43:23 $
 *    $Revision: 1.2.2.1 $
 *
 * Copyright (C) 2002-2007,
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
 


#include "DialogLoadVolume.h"
#include <wx/xrc/xmlres.h>
#include <wx/filedlg.h>
#include <wx/filename.h>

BEGIN_EVENT_TABLE( DialogLoadVolume, wxDialog )
	EVT_BUTTON			( wxID_OK,			 						DialogLoadVolume::OnOK )
	EVT_BUTTON			( XRCID( wxT( "ID_BUTTON_FILE" ) ),			DialogLoadVolume::OnButtonOpen )
	EVT_COMBOBOX		( XRCID( wxT( "ID_COMBO_FILENAME" ) ),		DialogLoadVolume::OnFileSelectionChanged )
END_EVENT_TABLE()


DialogLoadVolume::DialogLoadVolume( wxWindow* parent ) 
{
	wxXmlResource::Get()->LoadDialog( this, parent, wxT("ID_DIALOG_LOAD_VOLUME") );	
	m_checkNoResample = XRCCTRL( *this, "ID_CHECK_NO_RESAMPLE", wxCheckBox );
	m_btnOpen = XRCCTRL( *this, "ID_BUTTON_FILE", wxButton );
	m_comboFileName = XRCCTRL( *this, "ID_COMBO_FILENAME", wxComboBox );
	m_comboFileName->SetFocus();
}

DialogLoadVolume::~DialogLoadVolume()
{
}
	
wxString DialogLoadVolume::GetVolumeFileName()
{
	return m_comboFileName->GetValue().Trim( true ).Trim( false );
}

void DialogLoadVolume::OnOK( wxCommandEvent& event )
{
	if ( GetVolumeFileName().IsEmpty() )
	{
		wxMessageDialog dlg( this, "Volume file name can not be empty.", "Error", wxOK | wxICON_ERROR );
		dlg.ShowModal();
		return;
	}
		
	event.Skip();
}

bool DialogLoadVolume::IsToResample()
{
	return !m_checkNoResample->IsChecked();
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
	wxFileDialog dlg( this, _("Open volume file"), m_strLastDir, _(""), 
					  _T("Volume files (*.nii;*.nii.gz;*.img;*.mgz)|*.nii;*.nii.gz;*.img;*.mgz|All files (*.*)|*.*"), 
					  wxFD_OPEN );
	if ( dlg.ShowModal() == wxID_OK )
	{
		m_comboFileName->SetValue( dlg.GetPath() );
		m_comboFileName->SetInsertionPointEnd();
		m_strLastDir = wxFileName( dlg.GetPath() ).GetPath();
	}
}

void DialogLoadVolume::OnFileSelectionChanged( wxCommandEvent& event )
{
	m_comboFileName->SetInsertionPointEnd();
}
