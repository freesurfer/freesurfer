/**
 * @file  DialogLoadVolume.h
 * @brief Preferences Dialog.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2008/10/08 19:14:35 $
 *    $Revision: 1.5 $
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
#ifndef DialogLoadVolume_h
#define DialogLoadVolume_h

#include <wx/wx.h>

class DialogLoadVolume : public wxDialog
{
public:
	DialogLoadVolume( wxWindow* parent, bool bEnableResample = true );
	virtual ~DialogLoadVolume();
	
	wxString GetVolumeFileName();
	
	wxString GetRegFileName();
	
	bool IsToResample();
	
	void OnOK( wxCommandEvent& event );
	
	void SetLastDir( const wxString& dir )
		{ m_strLastDir = dir; }
	
	void SetRecentFiles( const wxArrayString& list );
	
protected:
	void OnButtonOpen( wxCommandEvent& event );
	void OnFileSelectionChanged( wxCommandEvent& event );
	void OnButtonRegFile( wxCommandEvent& event );
	void OnCheckApplyReg( wxCommandEvent& event );
	
	wxButton*		m_btnOpen;
	wxComboBox*		m_comboFileName;
	wxCheckBox*		m_checkNoResample;
	wxCheckBox*		m_checkApplyReg;
	wxTextCtrl*		m_textRegFile;
	wxButton*		m_btnRegFile;
	
	wxString		m_strLastDir;
	
	DECLARE_EVENT_TABLE()
};

#endif 

