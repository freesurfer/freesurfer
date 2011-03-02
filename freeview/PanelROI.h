/**
 * @file  PanelROI.h
 * @brief Layer control panel.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:02 $
 *    $Revision: 1.8 $
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
#ifndef PanelROI_h
#define PanelROI_h

#include <wx/wx.h>
#include "Listener.h"
#include "Broadcaster.h"


class wxAuiNotebook;
class wxListBox;
class wxCheckListBox;
class wxColourPickerCtrl;
class wxColourPickerEvent;
class Layer;

class PanelROI : public wxPanel, public Listener, public Broadcaster
{
public:
  PanelROI(wxWindow* parent);
  virtual ~PanelROI();

  void UpdateUI( bool bForce = false );

protected:
  void OnInternalIdle();

private:
  void OnSliderOpacity( wxScrollEvent& event );
  void OnLayerSelectionChanged( wxCommandEvent& event );
  void OnLayerVisibilityChanged( wxCommandEvent& event );

  void OnButtonNew( wxCommandEvent& event );
  void OnButtonLoad( wxCommandEvent& event );
  void OnButtonSave( wxCommandEvent& event );
  void OnButtonMoveUp( wxCommandEvent& event );
  void OnButtonMoveDown( wxCommandEvent& event );
  void OnButtonDelete( wxCommandEvent& event );
  void OnColorChanged( wxColourPickerEvent& event );
  void OnROICloseUpdateUI( wxUpdateUIEvent& event );
  void OnMoveUpUpdateUI( wxUpdateUIEvent& event );
  void OnMoveDownUpdateUI( wxUpdateUIEvent& event );

  void DoUpdateUI();

  void UpdateLayerList( Layer* layer );

  virtual void DoListenToMessage( std::string const iMsg, void* iData, void* sender );

  wxCheckListBox* m_listBoxLayers;
  wxButton*  m_btnMoveUp;
  wxButton*  m_btnMoveDown;
  wxButton*  m_btnNew;
  wxButton*  m_btnLoad;
  wxButton*  m_btnSave;
  wxButton*  m_btnDelete;
  wxSlider*  m_sliderOpacity;
  wxColourPickerCtrl*  m_colorPicker;
  wxTextCtrl*  m_textFileName;

  bool   m_bUINeedUpdate;

  DECLARE_EVENT_TABLE()
};

#endif

