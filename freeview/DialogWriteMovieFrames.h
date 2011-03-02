/**
 * @file  DialogWriteMovieFrames.h
 * @brief Dialog to save a volume.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 22:00:36 $
 *    $Revision: 1.4 $
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
#ifndef DialogWriteMovieFrames_h
#define DialogWriteMovieFrames_h

#include <wx/wx.h>

class wxComboBox;
class wxCheckBox;
class wxButton;

class DialogWriteMovieFrames : public wxDialog
{
public:
  DialogWriteMovieFrames( wxWindow* parent );
  virtual ~DialogWriteMovieFrames();

  wxString GetOutputLocation();
  
  wxString GetOutputExtension();
  
  double GetAngleStep();
  
  void UpdateUI();
  
protected:
  void OnWrite  ( wxCommandEvent& event );
  void OnAbort  ( wxCommandEvent& event );
  void OnClose  ( wxCommandEvent& event );
  void OnOpen   ( wxCommandEvent& event );
  
  void OnShow   ( wxShowEvent& event );
  

  wxTextCtrl*   m_textOutputLocation;
  wxChoice*     m_choiceOutputExtension;
  wxTextCtrl*   m_textAngleStep;
  wxButton*     m_buttonAbort;
  wxButton*     m_buttonWrite;
  wxButton*     m_buttonClose;
  
  DECLARE_EVENT_TABLE()
};

#endif

