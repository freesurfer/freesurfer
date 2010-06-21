/**
 * @file  DialogCropVolume.h
 * @brief Dialog to crop a volume.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2010/06/21 18:37:50 $
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
#ifndef DialogCropVolume_h
#define DialogCropVolume_h

#include <wx/wx.h>
#include "Listener.h"

class wxButton;
class wxSpinCtrl;
class LayerMRI;
class wxSpinEvent;

class DialogCropVolume : public wxDialog, public Listener
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
  
  virtual void DoListenToMessage ( std::string const iMessage, void* iData, void* sender );
  
  wxSpinCtrl*   m_spinRange[6];
  
  LayerMRI*     m_mri;
  
  DECLARE_EVENT_TABLE()
};

#endif

