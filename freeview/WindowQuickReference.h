/**
 * @file  WindowQuickReference.h
 * @brief Main window.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2008/06/04 20:43:25 $
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
 
#ifndef WindowQuickReference_h
#define WindowQuickReference_h

#include <wx/frame.h>

class wxHtmlWindow;

class WindowQuickReference : public wxFrame
{
public:
    WindowQuickReference( wxWindow* parent );
	virtual ~WindowQuickReference() {}
	
	void OnClose( wxCloseEvent& event );
    
protected:
	wxHtmlWindow*	m_wndHtml;
	
    // any class wishing to process wxWindows events must use this macro
	DECLARE_EVENT_TABLE()
};

#endif 


