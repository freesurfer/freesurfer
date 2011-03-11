/**
 * @file  PanelSceneSetting.h
 * @brief Layer control panel.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: krish $
 *    $Date: 2011/03/11 23:27:40 $
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
#ifndef PanelSceneSetting_h
#define PanelSceneSetting_h

#include <wx/wx.h>
#include "Listener.h"
#include "Broadcaster.h"


class wxAuiNotebook;
class wxListBox;
class wxCheckListBox;
class wxColourPickerCtrl;
class wxColourPickerEvent;
class Layer;

class PanelSceneSetting : public wxPanel, public Listener, public Broadcaster
{
public:
  PanelSceneSetting(wxWindow* parent);
  virtual ~PanelSceneSetting();

  void UpdateUI( bool bForce = false );

protected:
  void OnInternalIdle();

private:
  void OnSliderX( wxScrollEvent& event );
  void OnSliderY( wxScrollEvent& event );
  void OnSliderZ( wxScrollEvent& event );

  void OnCheckX( wxCommandEvent& event );
  void OnCheckY( wxCommandEvent& event );
  void OnCheckZ( wxCommandEvent& event );

  void DoUpdateUI();

  void UpdateLayerList( Layer* layer );

  virtual void DoListenToMessage( std::string const iMsg, void* iData, void* sender );

  wxCheckBox*  m_checkX;
  wxCheckBox*  m_checkY;
  wxCheckBox*  m_checkZ;

  bool   m_bUINeedUpdate;

  DECLARE_EVENT_TABLE()
};

#endif

