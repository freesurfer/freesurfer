/**
 * @file  DialogLoadDTI.h
 * @brief Preferences Dialog.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2009/04/29 22:53:49 $
 *    $Revision: 1.2.2.3 $
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
#ifndef DialogLoadDTI_h
#define DialogLoadDTI_h

#include <wx/wx.h>

class wxTextCtrl;
class wxCheckBox;
class wxComboBox;

class DialogLoadDTI : public wxDialog
{
public:
  DialogLoadDTI( wxWindow* parent );
  virtual ~DialogLoadDTI();

  wxString GetVectorFileName();
  wxString GetFAFileName();
  wxString GetRegFileName();

  bool IsToResample();

  void OnOK( wxCommandEvent& event );

  void SetLastDir( const wxString& dir )
  {
    m_strLastDir = dir;
  }

  void Initialize( bool bResample, bool bEnableCheckBox );

  void SetRecentFiles( const wxArrayString& list );

protected:
  void OnButtonVector( wxCommandEvent& event );
  void OnButtonFA( wxCommandEvent& event );
  void OnComboFASelectionChanged( wxCommandEvent& event );
  void OnButtonReg( wxCommandEvent& event );
  void OnCheckApplyReg( wxCommandEvent& event );

  wxButton*  m_btnVector;
  wxButton*  m_btnFA;
  wxButton*  m_btnReg;
  wxTextCtrl*  m_textVector;
  wxTextCtrl*  m_textReg;
  wxComboBox*  m_comboFA;
  wxCheckBox*  m_checkResample;
  wxCheckBox*  m_checkReg;

  wxString  m_strLastDir;

  DECLARE_EVENT_TABLE()
};

#endif

