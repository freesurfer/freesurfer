/**
 * @file  DialogRepositionSurface.h
 * @brief Dialog to reposition surface vertex.
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
#ifndef DialogRepositionSurface_h
#define DialogRepositionSurface_h

#include <wx/wx.h>
#include "Listener.h"

class wxTextCtrl;
class wxButton;
class wxCheckBox;
class wxChoice;
class LayerSurface;

class DialogRepositionSurface : public wxDialog, public Listener
{
public:
  DialogRepositionSurface( wxWindow* parent );
  virtual ~DialogRepositionSurface();

  void OnApply  ( wxCommandEvent& event );
  void OnUndo   ( wxCommandEvent& event );
  void OnSave   ( wxCommandEvent& event );
  void OnSaveAs ( wxCommandEvent& event );
  void OnClose  ( wxCommandEvent& event )
  {
    Close();
  }
  
  void OnChoiceTarget( wxCommandEvent& event );
  
  int GetVertex();
  int GetNeighborSize();
  
  double GetIntensity();
  double GetSigma();
  void GetCoordinate( double* pos );
  
  void UpdateUI();
  
protected:
  bool ValidateAll();
  void DoListenToMessage ( std::string const iMsg, void* iData, void* sender );

  wxTextCtrl*   m_textVertex;
  wxTextCtrl*   m_textTarget;
  wxTextCtrl*   m_textSize;
  wxTextCtrl*   m_textSigma;
  wxCheckBox*   m_checkUpdateVertex;
  wxChoice*     m_choiceTarget;
  wxButton*     m_btnSave;
  wxButton*     m_btnSaveAs;
  wxButton*     m_btnUndo;
  
  DECLARE_EVENT_TABLE()
};

#endif

