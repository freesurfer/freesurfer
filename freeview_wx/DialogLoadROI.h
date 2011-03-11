/**
 * @file  DialogLoadROI.h
 * @brief Dialog to load ROI/label data.
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
#ifndef DialogLoadROI_h
#define DialogLoadROI_h

#include <wx/wx.h>

class wxButton;
class wxChoice;
class wxTextCtrl;

class DialogLoadROI : public wxDialog
{
public:
  DialogLoadROI( wxWindow* parent, bool bEnableResample = true );
  virtual ~DialogLoadROI();

  wxArrayString GetFileNames();

  void OnOK( wxCommandEvent& event );

  void SetLastDir( const wxString& dir )
  {
    m_strLastDir = dir;
  }

  wxString GetTemplate();

protected:
  void OnButtonOpen           ( wxCommandEvent& event );

  wxTextCtrl*   m_textFileName;
  wxChoice*     m_choiceTemplate;

  wxString  m_strLastDir;

  DECLARE_EVENT_TABLE()
};

#endif

