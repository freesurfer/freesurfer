/**
 * @file  DialogCropVolume.h
 * @brief Dialog to crop a volume.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:01 $
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
#ifndef DialogCropVolume_h
#define DialogCropVolume_h

#include <wx/wx.h>
#include "Listener.h"

class wxButton;
class wxSpinCtrl;
class LayerMRI;
class wxSpinEvent;

class DialogCropVolume : public wxFrame, public Listener
{
public:
  DialogCropVolume( wxWindow* parent, LayerMRI* layer = NULL );
  virtual ~DialogCropVolume();

  void SetVolume( LayerMRI* mri );
  
protected:
  void OnButtonReset    ( wxCommandEvent& event );
  void OnButtonClose    ( wxCommandEvent& event );
  void OnButtonApply    ( wxCommandEvent& event );
  void OnButtonSaveAs   ( wxCommandEvent& event );
  void OnSpinBound      ( wxSpinEvent& event );
  void OnSpinBoundText  ( wxCommandEvent& event );
  void OnShow           ( wxShowEvent& event );
  void OnClose          ( wxCloseEvent& event );
  
  virtual void DoListenToMessage ( std::string const iMessage, void* iData, void* sender );
  
  wxSpinCtrl*   m_spinRange[6];
  
  LayerMRI*     m_mri;
  
  DECLARE_EVENT_TABLE()
};

#endif

