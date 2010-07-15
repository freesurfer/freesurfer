/**
 * @file  DialogTransformVolume.h
 * @brief Dialog to transform volume.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2010/07/15 19:51:47 $
 *    $Revision: 1.2 $
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
#ifndef DialogTransformVolume_h
#define DialogTransformVolume_h

#include <wx/wx.h>
#include "Listener.h"

class wxNotebook;
class wxNotebookEvent;

class DialogTransformVolume : public wxDialog, public Listener
{
public:
  DialogTransformVolume( wxWindow* parent );
  virtual ~DialogTransformVolume();

  void OnApply  ( wxCommandEvent& event );
  void OnRestore( wxCommandEvent& event );
  void OnSaveAs ( wxCommandEvent& event );
  void OnSaveReg( wxCommandEvent& event );

  bool GetRotation( int nID_in, int& plane_out, double& angle_out );

  void UpdateUI( int scope = 2 );
  
protected:
  void DoListenToMessage ( std::string const iMessage, void* iData, void* sender );
  
  void OnCheck1( wxCommandEvent& event );
  void OnCheck2( wxCommandEvent& event );
  void OnCheck3( wxCommandEvent& event );
  
  void OnTextTranslateX   ( wxCommandEvent& event );
  void OnTextTranslateY   ( wxCommandEvent& event );
  void OnTextTranslateZ   ( wxCommandEvent& event );
  void OnScrollTranslateX ( wxScrollEvent& event );
  void OnScrollTranslateY ( wxScrollEvent& event );
  void OnScrollTranslateZ ( wxScrollEvent& event );
  
  void OnPageChanged      ( wxNotebookEvent& event );
  
  void RespondTextTranslate   ( int n );
  void RespondScrollTranslate ( int n );
  
  void OnTextScaleX       ( wxCommandEvent& event );
  void OnTextScaleY       ( wxCommandEvent& event );
  void OnTextScaleZ       ( wxCommandEvent& event );
  void OnScrollScaleX     ( wxScrollEvent& event );
  void OnScrollScaleY     ( wxScrollEvent& event );
  void OnScrollScaleZ     ( wxScrollEvent& event );
  void RespondTextScale   ( int n );
  void RespondScrollScale ( int n );

  void DoRotate();

  wxCheckBox*   m_check[3];
  wxChoice*     m_choice[3];
  wxTextCtrl*   m_textAngle[3];
  wxButton*     m_btnRestoreOriginal;
  wxButton*     m_btnSaveAs;
  wxButton*     m_btnSaveReg;
  
  wxRadioButton*  m_radioAroundCenter;
  wxRadioButton*  m_radioAroundCursor;
  wxRadioButton*  m_radioNearest;
  wxRadioButton*  m_radioTrilinear;
  wxRadioButton*  m_radioSinc; 
  wxRadioButton*  m_radioActiveVolume;
  wxRadioButton*  m_radioAllVolumes;
  
  wxScrollBar*    m_scrollTranslate[3];
  wxTextCtrl*     m_textTranslate[3];
  wxScrollBar*    m_scrollScale[3];
  wxTextCtrl*     m_textScale[3];
  
  wxNotebook*     m_notebook;

  DECLARE_EVENT_TABLE()
};

#endif

