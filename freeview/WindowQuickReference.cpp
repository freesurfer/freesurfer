/**
 * @file  WindowQuickReference.h
 * @brief Main window.
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
 
#include "WindowQuickReference.h"
#include <wx/xrc/xmlres.h>
#include <wx/html/htmlwin.h>
#include <wx/config.h>
#include <wx/fs_mem.h>

BEGIN_EVENT_TABLE( WindowQuickReference, wxFrame )
	EVT_CLOSE		( WindowQuickReference::OnClose )
END_EVENT_TABLE()

#include "res/QuickRef.h"
WindowQuickReference::WindowQuickReference( wxWindow* parent )
{
	wxXmlResource::Get()->LoadFrame( this, parent, wxT("ID_FRAME_QUICK_REFERENCE") );
	m_wndHtml = XRCCTRL( *this, "ID_HTML_WINDOW", wxHtmlWindow );

	wxMemoryFSHandler::AddFileWithMimeType( "QuickRef.html", QuickRef_binary, QuickRef_binary_LEN, _T("text/html") );
	m_wndHtml->LoadPage( wxT("memory:QuickRef.html") );
//	m_wndHtml->LoadPage( "/homes/5/rpwang/freesurfer/dev/freeview/res/QuickRef.html" );
	wxConfigBase* config = wxConfigBase::Get();
	if ( config )
	{
		int w = config->Read( _T("/QuickRefWindow/Width"), 380L );
		int h = config->Read( _T("/QuickRefWindow/Height"), 520L );
		SetSize( w, h );
	}
	wxMemoryFSHandler::RemoveFile( "QuickRef.html" );
}

void WindowQuickReference::OnClose( wxCloseEvent& event )
{
	wxConfigBase* config = wxConfigBase::Get();
	if (config && !IsIconized())
	{
		int x, x2, y, y2, w, h;
		GetParent()->GetPosition( &x2, &y2 );
		GetPosition( &x, &y );
		GetSize( &w, &h );
		config->Write( _T("/QuickRefWindow/PosX"), (long) x - x2 );
		config->Write( _T("/QuickRefWindow/PosY"), (long) y - y2 );
		config->Write( _T("/QuickRefWindow/Width"), (long) w );
		config->Write( _T("/QuickRefWindow/Height"), (long) h );
	}
		
	Hide();
}
