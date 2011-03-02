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
 *    $Revision: 1.5 $
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

#ifndef WindowQuickReference_h
#define WindowQuickReference_h

#include <wx/frame.h>

class wxHtmlWindow;

class WindowQuickReference : public wxFrame
{
public:
  WindowQuickReference( wxWindow* parent );
  virtual ~WindowQuickReference()
  {}

  void OnClose( wxCloseEvent& event );

protected:
  wxHtmlWindow* m_wndHtml;

  // any class wishing to process wxWindows events must use this macro
  DECLARE_EVENT_TABLE()
};

#endif


