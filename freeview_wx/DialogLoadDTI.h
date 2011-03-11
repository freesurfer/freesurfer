/**
 * @file  DialogLoadDTI.h
 * @brief Dialog to load DTI data.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: krish $
 *    $Date: 2011/03/11 23:27:36 $
 *    $Revision: 1.1 $
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

