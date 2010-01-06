/**
 * @file  DialogPreferences.h
 * @brief Preferences dialog.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2010/01/06 02:07:19 $
 *    $Revision: 1.11 $
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



#include "DialogPreferences.h"
#include "MainWindow.h"
#include <wx/wx.h>
#include <wx/clrpicker.h>
#include <wx/xrc/xmlres.h>
#include <wx/spinctrl.h>
#include "stdlib.h"
#include "stdio.h"

BEGIN_EVENT_TABLE( DialogPreferences, wxDialog )
EVT_BUTTON     ( wxID_OK,         DialogPreferences::OnOK )
END_EVENT_TABLE()


DialogPreferences::DialogPreferences( wxWindow* parent )
{
  wxXmlResource::Get()->LoadDialog( this, parent, wxT("ID_DIALOG_PREFERENCES") );
  m_colorPickerBackground   = XRCCTRL( *this, "ID_COLORPICKER_BACKGROUND", wxColourPickerCtrl );
  m_colorPickerBackground->SetFocus();
  m_colorPickerCursor       = XRCCTRL( *this, "ID_COLORPICKER_CURSOR", wxColourPickerCtrl );
  m_choiceCursorStyle       = XRCCTRL( *this, "ID_CHOICE_CURSOR_STYLE", wxChoice );
  m_checkSaveCopy           = XRCCTRL( *this, "ID_CHECK_SAVE_COPY", wxCheckBox );
  m_checkSyncZoomFactor     = XRCCTRL( *this, "ID_CHECK_SYNC_ZOOM", wxCheckBox );

  m_checkHideCursor   = XRCCTRL( *this, "ID_CHECK_HIDE_CURSOR", wxCheckBox );
  m_checkHideCoords   = XRCCTRL( *this, "ID_CHECK_HIDE_COORDS", wxCheckBox );
  m_checkAntiAliasing = XRCCTRL( *this, "ID_CHECK_ANTIALIASING", wxCheckBox );  
  m_spinMagnification = XRCCTRL( *this, "ID_SPIN_MAGNIFICATION", wxSpinCtrl );
}

DialogPreferences::~DialogPreferences()
{}

void DialogPreferences::OnOK( wxCommandEvent& event )
{
  event.Skip();
}

void DialogPreferences::SetGeneralSettings( const SettingsGeneral& s )
{
  m_colorPickerBackground->SetColour( s.BackgroundColor );
  m_colorPickerCursor->SetColour( s.CursorColor );
  m_choiceCursorStyle->SetSelection( s.CursorStyle );
  m_checkSaveCopy->SetValue( s.SaveCopy );
}

SettingsGeneral DialogPreferences::GetGeneralSettings()
{
  SettingsGeneral s;
  s.BackgroundColor = m_colorPickerBackground->GetColour();
  s.CursorColor = m_colorPickerCursor->GetColour();
  s.CursorStyle = m_choiceCursorStyle->GetSelection();
  s.SaveCopy    = m_checkSaveCopy->GetValue();

  return s;
}


void DialogPreferences::Set2DSettings( const Settings2D& s )
{
  m_checkSyncZoomFactor->SetValue( s.SyncZoomFactor );
}

Settings2D DialogPreferences::Get2DSettings()
{
  Settings2D s;
  s.SyncZoomFactor = m_checkSyncZoomFactor->IsChecked();

  return s;
}

void DialogPreferences::SetScreenshotSettings( const SettingsScreenshot& s )
{
  m_spinMagnification ->SetValue( s.Magnification );
  m_checkHideCursor   ->SetValue( s.HideCursor );
  m_checkHideCoords   ->SetValue( s.HideCoords );
  m_checkAntiAliasing ->SetValue( s.AntiAliasing );
}

SettingsScreenshot DialogPreferences::GetScreenshotSettings()
{
  SettingsScreenshot s;
  s.HideCursor = m_checkHideCursor->IsChecked();
  s.HideCoords = m_checkHideCoords->IsChecked();
  s.Magnification = m_spinMagnification->GetValue();
  s.AntiAliasing = m_checkAntiAliasing->IsChecked();

  return s;
}
