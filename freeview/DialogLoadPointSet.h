/**
 * @file  DialogLoadPointSet.h
 * @brief Dialog to load DTI data.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2009/12/08 22:21:21 $
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

