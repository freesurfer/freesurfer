/**
 * @file  PanelSurface.h
 * @brief Layer control panel.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2009/01/27 18:27:25 $
 *    $Revision: 1.9 $
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
#ifndef PanelSurface_h
#define PanelSurface_h

#include <wx/wx.h>
#include "Listener.h"
#include "Broadcaster.h"
#include <vector>

class wxAuiNotebook;
class wxListBox;
class wxCheckListBox;
class wxColourPickerCtrl;
class wxColourPickerEvent;
class wxSpinCtrl;
class wxSpinEvent;
class wxTextCtrl;
class Layer;

class PanelSurface : public wxPanel, public Listener, public Broadcaster
{
public:
  PanelSurface(wxWindow* parent);
  virtual ~PanelSurface();

  void UpdateUI( bool bForce = false );

protected:
  void OnInternalIdle();

private:
  void OnSliderOpacity   ( wxScrollEvent& event );
  void OnLayerSelectionChanged ( wxCommandEvent& event );
  void OnLayerVisibilityChanged ( wxCommandEvent& event );

  void OnButtonLoad   ( wxCommandEvent& event );
  void OnButtonDelete   ( wxCommandEvent& event );
  void OnSurfaceCloseUpdateUI ( wxUpdateUIEvent& event );
  void OnColorChanged   ( wxColourPickerEvent& event );
  void OnEdgeColorChanged  ( wxColourPickerEvent& event );
  void OnSpinEdgeThickness ( wxSpinEvent& event );
  void OnButtonSurfaceMain ( wxCommandEvent& event );
  void OnButtonSurfaceInflated( wxCommandEvent& event );
  void OnButtonSurfaceWhite ( wxCommandEvent& event );
  void OnButtonSurfacePial ( wxCommandEvent& event );
  void OnButtonSurfaceOriginal( wxCommandEvent& event );
  void OnChoiceVector   ( wxCommandEvent& event );
  void OnVectorColorChanged ( wxColourPickerEvent& event );
  void OnSpinVectorPointSize ( wxSpinEvent& event );

  void OnChoiceCurvatureMap ( wxCommandEvent& evet );
  void OnSliderMidPoint  ( wxScrollEvent& event );
  void OnSliderSlope   ( wxScrollEvent& event );
  void OnSliderMidPointChanging  ( wxScrollEvent& event );
  void OnSliderSlopeChanging   ( wxScrollEvent& event );
  void OnTextMidPoint   ( wxCommandEvent& event );
  void OnTextSlope   ( wxCommandEvent& event );

  void DoUpdateUI();

  void UpdateTextValue( wxTextCtrl* ctrl, double dvalue );

  void UpdateLayerList( Layer* layer );

  virtual void DoListenToMessage( std::string const iMsg, void* iData );

  wxCheckListBox* m_listBoxLayers;
  wxButton*  m_btnNew;
  wxButton*  m_btnLoad;
  wxButton*  m_btnSave;
  wxButton*  m_btnDelete;
  wxButton*  m_btnSurfaceMain;
  wxButton*  m_btnSurfaceInflated;
  wxButton*  m_btnSurfaceWhite;
  wxButton*  m_btnSurfacePial;
  wxButton*  m_btnSurfaceOriginal;
  wxSlider*  m_sliderOpacity;
  wxColourPickerCtrl*  m_colorPicker;
  wxColourPickerCtrl*  m_colorPickerEdge;
  wxTextCtrl*  m_textFileName;
  wxSpinCtrl*  m_spinEdgeThickness;
  wxChoice*  m_choiceVector;
  wxColourPickerCtrl*  m_colorPickerVector;
  wxSpinCtrl*  m_spinVectorPointSize;

  wxChoice*  m_choiceCurvatureMap;
  wxSlider*  m_sliderMidPoint;
  wxTextCtrl*  m_textMidPoint;
  wxSlider*  m_sliderSlope;
  wxTextCtrl*   m_textSlope;

  std::vector<wxWindow*> m_widgetsMidPoint;
  std::vector<wxWindow*> m_widgetsSlope;
  std::vector<wxWindow*> m_widgetsVector;

  bool   m_bUINeedUpdate;

  DECLARE_EVENT_TABLE()
};

#endif

