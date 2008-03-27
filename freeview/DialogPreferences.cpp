/**
 * @file  DialogPreferences.h
 * @brief Preferences dialog.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2008/03/27 18:12:14 $
 *    $Revision: 1.1 $
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
	m_colorPicker = XRCCTRL( *this, "ID_COLORPICKER", wxColourPickerCtrl );
	m_colorPicker->SetFocus();
}

DialogPreferences::~DialogPreferences()
{
}

wxColour DialogPreferences::GetBackgroundColor() const
{
	return m_colorPicker->GetColour();
}

void DialogPreferences::SetBackgroundColor( const wxColour& color )
{
	m_colorPicker->SetColour( color );
}

void DialogPreferences::OnOK( wxCommandEvent& event )
{
	
	event.Skip();
}
