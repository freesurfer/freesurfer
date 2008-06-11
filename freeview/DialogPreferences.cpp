/**
 * @file  DialogPreferences.h
 * @brief Preferences dialog.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2008/06/11 21:30:18 $
 *    $Revision: 1.3 $
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
