/**
 * @file  PanelSurface.h
 * @brief Layer control panel.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2009/10/20 21:41:40 $
 *    $Revision: 1.16 $
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
  void OnSurfaceClose           ( wxCommandEvent& event );
  void OnSurfaceCloseUpdateUI   ( wxUpdateUIEvent& event );
  
  void OnSliderOpacityChanging  ( wxScrollEvent& event );
  void OnSliderOpacity          ( wxScrollEvent& event );
  void OnTextOpacity            ( wxCommandEvent& event );
  void OnLayerSelectionChanged  ( wxCommandEvent& event );
  void OnLayerVisibilityChanged ( wxCommandEvent& event );
  
  void OnColorChanged           ( wxColourPickerEvent& event );
  void OnEdgeColorChanged       ( wxColourPickerEvent& event );
  void OnSpinEdgeThickness      ( wxSpinEvent& event );

  void OnChoiceVector           ( wxCommandEvent& event );
  void OnVectorColorChanged     ( wxColourPickerEvent& event );
  void OnSpinVectorPointSize    ( wxSpinEvent& event );
  
  void OnChoiceOverlay          ( wxCommandEvent& event );
  void OnButtonConfigureOverlay ( wxCommandEvent& event );          

  void OnChoiceCurvatureMap     ( wxCommandEvent& evet );
  void OnSliderMidPoint         ( wxScrollEvent& event );
  void OnSliderSlope            ( wxScrollEvent& event );
  void OnSliderMidPointChanging ( wxScrollEvent& event );
  void OnSliderSlopeChanging    ( wxScrollEvent& event );
  void OnTextMidPoint           ( wxCommandEvent& event );
  void OnTextSlope              ( wxCommandEvent& event );
  void OnChoiceRenderMode       ( wxCommandEvent& event );
  void OnChoiceMeshColorMap     ( wxCommandEvent& event );
  
  void OnChoiceAnnotation       ( wxCommandEvent& event );
  
  void OnCheckShowVertices      ( wxCommandEvent& event );
  void OnColorVertex            ( wxColourPickerEvent& event );
  void OnSpinVertexPointSize    ( wxSpinEvent& event );
  
  void OnTextPosition           ( wxCommandEvent& event );

  void DoUpdateUI();

  void UpdateTextValue( wxTextCtrl* ctrl, double dvalue );

  void UpdateLayerList( Layer* layer );

  virtual void DoListenToMessage( std::string const iMsg, void* iData, void* sender );

  wxCheckListBox* m_listBoxLayers;
  wxButton*     m_btnConfigureOverlay;
  wxSlider*     m_sliderOpacity;
  wxTextCtrl*   m_textOpacity;
  wxColourPickerCtrl*  m_colorPicker;
  wxColourPickerCtrl*  m_colorPickerEdge;
  wxTextCtrl*   m_textFileName;
  wxSpinCtrl*   m_spinEdgeThickness;
  wxChoice*     m_choiceVector;
  wxColourPickerCtrl*  m_colorPickerVector;
  wxSpinCtrl*   m_spinVectorPointSize;

  wxChoice*     m_choiceCurvatureMap;
  wxSlider*     m_sliderMidPoint;
  wxTextCtrl*   m_textMidPoint;
  wxSlider*     m_sliderSlope;
  wxTextCtrl*   m_textSlope;
  
  wxChoice*     m_choiceOverlay;
  wxButton*     m_btnOverlayConfiguration;
  
  wxChoice*     m_choiceAnnotation;
  
  wxChoice*     m_choiceRenderMode;
  wxChoice*     m_choiceMeshColorMap;
  
  wxCheckBox*   m_checkShowVertices;
  wxColourPickerCtrl* m_colorPickerVertex;
  wxSpinCtrl*   m_spinVertexPointSize;
  
  wxTextCtrl*   m_textPosition;

  std::vector<wxWindow*>  m_widgetsMidPoint;
  std::vector<wxWindow*>  m_widgetsSlope;
  std::vector<wxWindow*>  m_widgetsVector;
  std::vector<wxWindow*>  m_widgetsVertex;
  std::vector<wxWindow*>  m_widgetsMesh;

  bool   m_bUINeedUpdate;

  DECLARE_EVENT_TABLE()
};

#endif

