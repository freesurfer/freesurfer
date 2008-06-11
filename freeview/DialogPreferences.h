/**
 * @file  DialogPreferences.h
 * @brief Preferences Dialog.
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
#ifndef DialogPreferences_h
#define DialogPreferences_h

#include <wx/wx.h>


class wxColourPickerCtrl;

class DialogPreferences : public wxDialog
{
public:
	DialogPreferences(wxWindow* parent);
	virtual ~DialogPreferences();
	
	wxColour GetBackgroundColor() const;
	void SetBackgroundColor( const wxColour& color );
	
	wxColour GetCursorColor() const;
	void SetCursorColor( const wxColour& color );
	
	void OnOK( wxCommandEvent& event ); 
	
private:
	wxColourPickerCtrl*		m_colorPickerBackground;
	wxColourPickerCtrl*		m_colorPickerCursor;
	
	DECLARE_EVENT_TABLE()
};

#endif 

