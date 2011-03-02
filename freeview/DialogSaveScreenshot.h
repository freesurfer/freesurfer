/**
 * @file  DialogSaveScreenshot.h
 * @brief Dialog to save screenshot of the current view.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:01 $
 *    $Revision: 1.2 $
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

