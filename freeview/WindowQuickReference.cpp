/**
 * @file  WindowQuickReference.h
 * @brief Main window.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:03 $
 *    $Revision: 1.9 $
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

#include "WindowQuickReference.h"
#include <wx/xrc/xmlres.h>
#include <wx/html/htmlwin.h>
#include <wx/config.h>
#include <wx/fs_mem.h>

BEGIN_EVENT_TABLE( WindowQuickReference, wxFrame )
  EVT_CLOSE  ( WindowQuickReference::OnClose )
END_EVENT_TABLE()

#include "res/QuickRef.h"
WindowQuickReference::WindowQuickReference( wxWindow* parent )
{
  wxXmlResource::Get()->LoadFrame( this, parent, wxT("ID_FRAME_QUICK_REFERENCE") );
  m_wndHtml = XRCCTRL( *this, "ID_HTML_WINDOW", wxHtmlWindow );

  wxMemoryFSHandler::AddFileWithMimeType( _("QuickRef.html"), QuickRef_binary, QuickRef_binary_LEN, _("text/html") );
  m_wndHtml->LoadPage( _("memory:QuickRef.html") );
  wxConfigBase* config = wxConfigBase::Get();
  if ( config )
  {
    int w = config->Read( _("/QuickRefWindow/Width"), 380L );
    int h = config->Read( _("/QuickRefWindow/Height"), 520L );
    SetSize( w, h );
  }
  wxMemoryFSHandler::RemoveFile( _("QuickRef.html") );
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
    config->Write( _("/QuickRefWindow/PosX"),   (long) x - x2 );
    config->Write( _("/QuickRefWindow/PosY"),   (long) y - y2 );
    config->Write( _("/QuickRefWindow/Width"),  (long) w );
    config->Write( _("/QuickRefWindow/Height"), (long) h );
  }

  Hide();
}
