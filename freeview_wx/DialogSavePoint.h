/**
 * @file  DialogSavePoint.h
 * @brief Dialog to save a volume.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: krish $
 *    $Date: 2011/03/11 23:27:37 $
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
#ifndef DialogSavePoint_h
#define DialogSavePoint_h

#include <wx/wx.h>

class wxComboBox;
class wxCheckBox;

class DialogSavePoint : public wxDialog
{
public:
  DialogSavePoint( wxWindow* parent );
  virtual ~DialogSavePoint();

  wxString GetGroupName();
  wxString GetPointName();
  
  void SetGotoPoints( const wxArrayString& points );
  wxArrayString GetGotoPoints()
  {
    return m_points;
  }
  
protected:
  void OnSave   ( wxCommandEvent& event );
  void OnDelete ( wxCommandEvent& event );
  void OnClose  ( wxCommandEvent& event );
  void OnComboGroupName ( wxCommandEvent& event );
  
  void OnShow   ( wxShowEvent& event );
  
  void UpdateNames();
  bool ValidateInput();

  wxComboBox*   m_comboGroupName;
  wxComboBox*   m_comboPointName;

  wxArrayString m_points;
  
  DECLARE_EVENT_TABLE()
};

#endif

