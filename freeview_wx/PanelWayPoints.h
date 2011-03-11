/**
 * @file  PanelWayPoints.h
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
#ifndef PanelWayPoints_h
#define PanelWayPoints_h

#include <wx/wx.h>
#include "Listener.h"
#include "Broadcaster.h"
#include <vector>

class wxAuiNotebook;
class wxListBox;
class wxCheckListBox;
class wxCheckBox;
class wxColourPickerCtrl;
class wxColourPickerEvent;
class wxTextCtrl;
class wxSlider;
class wxChoice;
class Layer;

class PanelWayPoints : public wxPanel, public Listener, public Broadcaster
{
public:
  PanelWayPoints(wxWindow* parent);
  virtual ~PanelWayPoints();

  void UpdateUI( bool bForce = false );

protected:
  void OnInternalIdle();

private:
  void OnSliderOpacity   ( wxScrollEvent& event );
  void OnLayerSelectionChanged ( wxCommandEvent& event );
  void OnLayerVisibilityChanged ( wxCommandEvent& event );

  void OnWayPointsLoad   ( wxCommandEvent& event );
  void OnWayPointsClose   ( wxCommandEvent& event );
  void OnWayPointsCloseUpdateUI ( wxUpdateUIEvent& event );
  void OnColorChanged    ( wxColourPickerEvent& event );
  void OnSplineColorChanged  ( wxColourPickerEvent& event );
  void OnTextRadiusChanged  ( wxCommandEvent& event );
  void OnTextSplineRadius   ( wxCommandEvent& event );

  void OnChoiceColorMap   ( wxCommandEvent& event );
  void OnChoiceScalarMap   ( wxCommandEvent& event );
  void OnSliderHeatScaleMin  ( wxScrollEvent& event );
  void OnSliderHeatScaleMid  ( wxScrollEvent& event );
  void OnSliderHeatScaleMax  ( wxScrollEvent& event );
  void OnSliderHeatScaleOffset ( wxScrollEvent& event );
  void OnTextHeatScaleMin   ( wxCommandEvent& event );
  void OnTextHeatScaleMid   ( wxCommandEvent& event );
  void OnTextHeatScaleMax   ( wxCommandEvent& event );
  void OnTextHeatScaleOffset  ( wxCommandEvent& event );
  void OnCheckShowSpline      ( wxCommandEvent& event );
  void OnCheckSnapToVoxelCenter ( wxCommandEvent& event );

  void DoUpdateUI();

  void UpdateLayerList( Layer* layer );

  virtual void DoListenToMessage( std::string const iMsg, void* iData, void* sender );

  void UpdateTextValue( wxTextCtrl* ctrl, double dvalue );
  void ShowWidgets( std::vector<wxWindow*>& list, bool bShow );

  void LoadScalarValues();

  wxCheckListBox* m_listBoxLayers;
  wxSlider*  m_sliderOpacity;
  wxColourPickerCtrl*  m_colorPicker;
  wxColourPickerCtrl*  m_colorPickerSpline;
  wxTextCtrl* m_textFileName;
  wxTextCtrl* m_textRadius;
  wxTextCtrl* m_textSplineRadius;
  wxChoice*   m_choiceColorMap;
  wxChoice*   m_choiceScalarMap;
  wxSlider*   m_sliderHeatScaleMin;
  wxSlider*   m_sliderHeatScaleMid;
  wxSlider*   m_sliderHeatScaleMax;
  wxSlider*   m_sliderHeatScaleOffset;
  wxTextCtrl* m_textHeatScaleMin;
  wxTextCtrl* m_textHeatScaleMid;
  wxTextCtrl* m_textHeatScaleMax;
  wxTextCtrl* m_textHeatScaleOffset;
  wxCheckBox* m_checkShowSpline;
  wxCheckBox* m_checkSnapToVoxelCenter;

  std::vector<wxWindow*> m_widgetlistSolidColor;
  std::vector<wxWindow*> m_widgetlistHeatScale;
  std::vector<wxWindow*> m_widgetlistSpline;

  bool   m_bUINeedUpdate;

  DECLARE_EVENT_TABLE()
};

#endif

