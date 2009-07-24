/**
 * @file  DialogLoadVolume.h
 * @brief Dialog to load volume data.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2009/07/24 20:38:22 $
 *    $Revision: 1.12 $
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
#ifndef DialogLoadVolume_h
#define DialogLoadVolume_h

#include <wx/wx.h>

class wxButton;
class wxComboBox;
class wxTextCtrl;
class wxRadioButton;

class DialogLoadVolume : public wxDialog
{
public:
  DialogLoadVolume( wxWindow* parent, bool bEnableResample = true );
  virtual ~DialogLoadVolume();

  wxArrayString GetVolumeFileNames();

  wxString GetRegFileName();

  bool IsToResample();
  
  int GetSampleMethod();

  void OnOK( wxCommandEvent& event );

  void SetLastDir( const wxString& dir )
  {
    m_strLastDir = dir;
  }

  void SetRecentFiles( const wxArrayString& list );
  
  wxString GetColorMap();
  wxString GetLUT();

protected:
  void OnButtonOpen           ( wxCommandEvent& event );
  void OnFileSelectionChanged ( wxCommandEvent& event );
  void OnButtonRegFile        ( wxCommandEvent& event );
  void OnCheckApplyReg        ( wxCommandEvent& event );
  void OnChoiceColorMap       ( wxCommandEvent& event );
  void OnChoiceLUT            ( wxCommandEvent& event );
  
  void UpdateLUT();

  wxButton*     m_btnOpen;
  wxComboBox*   m_comboFileName;
  wxCheckBox*   m_checkResample;
  wxCheckBox*   m_checkApplyReg;
  wxTextCtrl*   m_textRegFile;
  wxButton*     m_btnRegFile;
  wxRadioButton*  m_radioNearest;
  wxRadioButton*  m_radioTrilinear;
  wxChoice*     m_choiceColorMap;
  wxChoice*     m_choiceLUT;
  wxStaticText* m_staticLUT;

  wxString  m_strLastDir;

  DECLARE_EVENT_TABLE()
};

#endif

