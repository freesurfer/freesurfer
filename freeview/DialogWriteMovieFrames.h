/**
 * @file  DialogWriteMovieFrames.h
 * @brief Dialog to save a volume.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2010/01/11 21:30:14 $
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
  
protected:
  void OnWrite  ( wxCommandEvent& event );
  void OnAbort  ( wxCommandEvent& event );
  void OnClose  ( wxCommandEvent& event );
  void OnOpen   ( wxCommandEvent& event );
  
  void OnShow   ( wxShowEvent& event );
  
  void UpdateUI ();

  wxTextCtrl*   m_textOutputLocation;
  wxChoice*     m_choiceOutputExtension;
  wxTextCtrl*   m_textAngleStep;
  wxButton*     m_buttonAbort;
  wxButton*     m_buttonWrite;
  wxButton*     m_buttonClose;
  
  DECLARE_EVENT_TABLE()
};

#endif

