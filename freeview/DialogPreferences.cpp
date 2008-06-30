/**
 * @file  DialogPreferences.h
 * @brief Preferences dialog.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2008/06/30 20:48:35 $
 *    $Revision: 1.4 $
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
 


#include "DialogPreferences.h"
#include "MainWindow.h"
#include <wx/wx.h>
#include <wx/clrpicker.h>
#include <wx/xrc/xmlres.h>
#include "stdlib.h"
#include "stdio.h"

BEGIN_EVENT_TABLE( DialogPreferences, wxDialog )
	EVT_BUTTON					( wxID_OK,			 					DialogPreferences::OnOK )
END_EVENT_TABLE()


DialogPreferences::DialogPreferences( wxWindow* parent ) 
{
	wxXmlResource::Get()->LoadDialog( this, parent, wxT("ID_DIALOG_PREFERENCES") );
	m_colorPickerBackground = XRCCTRL( *this, "ID_COLORPICKER_BACKGROUND", wxColourPickerCtrl );
	m_colorPickerBackground->SetFocus();
	m_colorPickerCursor = XRCCTRL( *this, "ID_COLORPICKER_CURSOR", wxColourPickerCtrl );
	m_checkSyncZoomFactor = XRCCTRL( *this, "ID_CHECK_SYNC_ZOOM", wxCheckBox );
}

DialogPreferences::~DialogPreferences()
{
}

wxColour DialogPreferences::GetBackgroundColor() const
{
	return m_colorPickerBackground->GetColour();
}

void DialogPreferences::SetBackgroundColor( const wxColour& color )
{
	m_colorPickerBackground->SetColour( color );
}


wxColour DialogPreferences::GetCursorColor() const
{
	return m_colorPickerCursor->GetColour();
}

void DialogPreferences::SetCursorColor( const wxColour& color )
{
	m_colorPickerCursor->SetColour( color );
}

void DialogPreferences::OnOK( wxCommandEvent& event )
{
	event.Skip();
}

void DialogPreferences::Set2DSettings( const Settings2D& s )
{
	m_checkSyncZoomFactor->SetValue( s.SyncZoomFactor );
}

Settings2D DialogPreferences::Get2DSettings()
{
	Settings2D s;
	s.SyncZoomFactor = m_checkSyncZoomFactor->IsChecked();
	
	return s;
}
