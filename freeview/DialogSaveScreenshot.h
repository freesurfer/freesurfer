/**
 * @file  DialogSaveScreenshot.h
 * @brief Dialog to save screenshot of the current view.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2009/12/04 21:57:12 $
 *    $Revision: 1.1 $
 *
 * Copyright (C) 2008-2009,
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
#ifndef DialogSaveScreenshot_h
#define DialogSaveScreenshot_h

#include <wx/wx.h>
#include "CommonDataStruct.h"

class wxTextCtrl;
class wxCheckBox;
class wxSpinCtrl;

class DialogSaveScreenshot : public wxDialog
{
public:
  DialogSaveScreenshot( wxWindow* parent );
  virtual ~DialogSaveScreenshot();

  wxString GetFileName();
  void SetFileName( const wxString& filename );

  void SetSettings( SettingsScreenshot s );
  
  SettingsScreenshot GetSettings();
  
  void OnOK( wxCommandEvent& event );

  void SetLastDir( const wxString& dir )
  {
    m_strLastDir = dir;
  }
  
  int GetFilterIndex()
  {
    return m_nScreenshotFilterIndex;
  }
  
  void SetFilterIndex( int nIndex )
  {
    m_nScreenshotFilterIndex = nIndex;
  }
  
protected:
  void OnButtonOpen( wxCommandEvent& event );

  wxButton*     m_btnOpen;
  wxTextCtrl*   m_textFilename;
  wxCheckBox*   m_checkAntiAliasing;
  wxCheckBox*   m_checkHideCursor;
  wxCheckBox*   m_checkHideCoords;
  wxCheckBox*   m_checkKeepWindow;
  wxSpinCtrl*   m_spinMagnification;
  
  int       m_nScreenshotFilterIndex;
  wxString  m_strLastDir;

  DECLARE_EVENT_TABLE()
};

#endif

