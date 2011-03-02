/**
 * @file  DialogPreferences.h
 * @brief Preferences Dialog.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 22:00:36 $
 *    $Revision: 1.12 $
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
#ifndef DialogPreferences_h
#define DialogPreferences_h

#include <wx/wx.h>


class wxColourPickerCtrl;
class wxCheckBox;
class wxSpinCtrl;
struct SettingsGeneral;
struct Settings2D;
struct SettingsScreenshot;

class DialogPreferences : public wxDialog
{
public:
  DialogPreferences(wxWindow* parent);
  virtual ~DialogPreferences();

  SettingsGeneral GetGeneralSettings();
  void SetGeneralSettings( const SettingsGeneral& s );

  Settings2D Get2DSettings();
  void Set2DSettings( const Settings2D& s );

  SettingsScreenshot GetScreenshotSettings();
  void SetScreenshotSettings( const SettingsScreenshot& s );

  void OnOK( wxCommandEvent& event );

private:
  wxColourPickerCtrl* m_colorPickerBackground;
  wxColourPickerCtrl* m_colorPickerCursor;
  wxChoice*           m_choiceCursorStyle;
  wxCheckBox*         m_checkSaveCopy;
  wxCheckBox*         m_checkSyncZoomFactor;

  wxCheckBox*         m_checkHideCursor;
  wxCheckBox*         m_checkHideCoords;
  wxCheckBox*         m_checkAntiAliasing;
  wxSpinCtrl*         m_spinMagnification;

  DECLARE_EVENT_TABLE()
};

#endif

