/**
 * @file  DialogLoadPointSet.h
 * @brief Dialog to load DTI data.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 22:00:36 $
 *    $Revision: 1.3 $
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
#ifndef DialogLoadPointSet_h
#define DialogLoadPointSet_h

#include <wx/wx.h>

class wxTextCtrl;
class wxCheckBox;
class wxRadioButton;

class DialogLoadPointSet : public wxDialog
{
public:
  DialogLoadPointSet( wxWindow* parent );
  virtual ~DialogLoadPointSet();

  wxString GetFileName();

  int GetPointSetType();

  void OnOK( wxCommandEvent& event );

  void SetLastDir( const wxString& dir )
  {
    m_strLastDir = dir;
  }

protected:
  void OnButtonOpen( wxCommandEvent& event );

  wxButton*       m_btnOpen;
  wxTextCtrl*     m_textFileName;
  wxRadioButton*  m_radioWayPoints;
  wxRadioButton*  m_radioControlPoints;

  wxString  m_strLastDir;

  DECLARE_EVENT_TABLE()
};

#endif

