/**
 * @file  DialogSavePointSetAs.h
 * @brief Dialog to save point set.
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
#ifndef DialogSavePointSetAs_h
#define DialogSavePointSetAs_h

#include <wx/wx.h>

class wxTextCtrl;
class wxCheckBox;
class wxRadioButton;

class DialogSavePointSetAs : public wxDialog
{
public:
  DialogSavePointSetAs( wxWindow* parent );
  virtual ~DialogSavePointSetAs();

  wxString GetFileName();

  void OnOK( wxCommandEvent& event );

  void SetLastDir( const wxString& dir )
  {
    m_strLastDir = dir;
  }
  
  void SetFileName( const wxString& filename );
  
  void SetPointSetType( int nType );
  
  int GetPointSetType();
  
protected:
  void OnButtonOpen( wxCommandEvent& event );

  wxButton*       m_btnOpen;
  wxTextCtrl*     m_textFilename;
  wxRadioButton*  m_radioWayPoints;
  wxRadioButton*  m_radioControlPoints;

  wxString  m_strLastDir;

  DECLARE_EVENT_TABLE()
};

#endif

