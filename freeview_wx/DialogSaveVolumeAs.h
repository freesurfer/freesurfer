/**
 * @file  DialogSaveVolumeAs.h
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
#ifndef DialogSaveVolumeAs_h
#define DialogSaveVolumeAs_h

#include <wx/wx.h>

class wxTextCtrl;
class wxCheckBox;

class DialogSaveVolumeAs : public wxDialog
{
public:
  DialogSaveVolumeAs( wxWindow* parent );
  virtual ~DialogSaveVolumeAs();

  wxString GetFileName();
  bool GetReorient();

  void OnOK( wxCommandEvent& event );

  void SetLastDir( const wxString& dir )
  {
    m_strLastDir = dir;
  }
  
  void SetFileName( const wxString& filename );
  
protected:
  void OnButtonOpen( wxCommandEvent& event );

  wxButton*  m_btnOpen;
  wxTextCtrl*  m_textFilename;
  wxCheckBox*  m_checkReorient;

  wxString  m_strLastDir;

  DECLARE_EVENT_TABLE()
};

#endif

