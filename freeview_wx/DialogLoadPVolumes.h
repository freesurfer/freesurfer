/**
 * @file  DialogLoadPVolumes.h
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
#ifndef DialogLoadPVolumes_h
#define DialogLoadPVolumes_h

#include <wx/wx.h>

class wxTextCtrl;
class wxCheckBox;
class wxChoice;

class DialogLoadPVolumes : public wxDialog
{
public:
  DialogLoadPVolumes( wxWindow* parent );
  virtual ~DialogLoadPVolumes();

  wxArrayString GetVolumeFileNames();
  
  wxString GetLUT();

  wxString GetFileNamePrefix();
  
  void SetLastDir( const wxString& dir )
  {
    m_strLastDir = dir;
  }

//  void Initialize( bool bResample, bool bEnableCheckBox );

protected:
  void OnOK             ( wxCommandEvent& event );
  void OnButtonOpen     ( wxCommandEvent& event );
  void OnChoiceLUT      ( wxCommandEvent& event );
  
  void UpdateLUT();

  wxButton*     m_btnOpen;
  wxTextCtrl*   m_textFileNames;
  wxTextCtrl*   m_textPrefix;
  wxChoice*     m_choiceLUT;

  wxString  m_strLastDir;

  DECLARE_EVENT_TABLE()
};

#endif

