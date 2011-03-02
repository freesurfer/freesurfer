/**
 * @file  DialogLoadVolume.h
 * @brief Dialog to load volume data.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:01 $
 *    $Revision: 1.13 $
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

