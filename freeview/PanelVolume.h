/**
 * @file  PanelVolume.h
 * @brief Layer control panel.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 22:00:36 $
 *    $Revision: 1.33 $
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
#ifndef PanelVolume_h
#define PanelVolume_h

#include <wx/wx.h>
#include "Listener.h"
#include "Broadcaster.h"
#include <vector>
#include "LUTDataHolder.h"

class wxAuiNotebook;
class wxListBox;
class wxCheckListBox;
class wxComboBox;
class wxChoice;
class wxTextCtrl;
class wxColorIndicator;
class wxColourPickerCtrl;
class wxColourPickerEvent;
class wxButton;
class Layer;
class LayerMRI;

class PanelVolume : public wxPanel, public Listener, public Broadcaster
{
public:
  PanelVolume(wxWindow* parent = NULL);
  virtual ~PanelVolume();

  void UpdateUI( bool bForce = false );

private:
  void InitWidgetsFromXRC(wxWindow* parent);
  void OnInternalIdle();

private:
  void OnSliderOpacityChanged   ( wxScrollEvent& event );
  void OnTextOpacityChanged     ( wxCommandEvent& event );
  void OnLayerSelectionChanged  ( wxCommandEvent& event );
  void OnLayerVisibilityChanged ( wxCommandEvent& event );
  void OnListDoubleClicked      ( wxCommandEvent& event );

  void OnButtonNew              ( wxCommandEvent& event );
  void OnButtonLoad             ( wxCommandEvent& event );
  void OnButtonSave             ( wxCommandEvent& event );
  void OnButtonMoveUp           ( wxCommandEvent& event );
  void OnMoveUpUpdateUI         ( wxUpdateUIEvent& event );
  void OnButtonMoveDown         ( wxCommandEvent& event );
  void OnMoveDownUpdateUI       ( wxUpdateUIEvent& event );
  void OnButtonDelete           ( wxCommandEvent& event );
  void OnVolumeCloseUpdateUI    ( wxUpdateUIEvent& event );
  void OnVolumeLock             ( wxCommandEvent& event );
  void OnVolumeLockUpdateUI     ( wxUpdateUIEvent& event );
  void OnVolumeCopySetting      ( wxCommandEvent& event );
  void OnVolumeCopySettingUpdateUI    ( wxUpdateUIEvent& event );
  void OnVolumePasteSetting     ( wxCommandEvent& event );
  void OnVolumePasteSettingUpdateUI   ( wxUpdateUIEvent& event );
  void OnVolumePasteSettingAll  ( wxCommandEvent& event );
  void OnVolumePasteSettingAllUpdateUI( wxUpdateUIEvent& event );

  void OnCheckClearBackground   ( wxCommandEvent& event );
  void OnChoiceColorMap         ( wxCommandEvent& event );
  void OnChoiceLUT              ( wxCommandEvent& event );
  void OnColorSelectionChanged  ( wxCommandEvent& event );
  void OnTextFillValueChanged   ( wxCommandEvent& event );
  void OnTextWindowChanged      ( wxCommandEvent& event );
  void OnTextLevelChanged       ( wxCommandEvent& event );
  void OnSliderWindowChanged    ( wxScrollEvent& event );
  void OnSliderLevelChanged     ( wxScrollEvent& event );
  void OnTextColorMapMinChanged       ( wxCommandEvent& event );
  void OnTextColorMapMaxChanged       ( wxCommandEvent& event );
  void OnTextHeatScaleMidChanged      ( wxCommandEvent& event );
  void OnTextHeatScaleOffsetChanged   ( wxCommandEvent& event );
  void OnSliderColorMapMinChanged     ( wxScrollEvent& event );
  void OnSliderColorMapMaxChanged     ( wxScrollEvent& event );
  void OnSliderHeatScaleMidChanged    ( wxScrollEvent& event );
  void OnSliderHeatScaleOffsetChanged ( wxScrollEvent& event );
  void OnCheckUpsample                ( wxCommandEvent& event );
  void OnCheckSmooth                  ( wxCommandEvent& event );
  void OnChoiceDirectionCode          ( wxCommandEvent& event );
  void OnSliderFrameChanged           ( wxScrollEvent& event );
  void OnTextFrameChanged             ( wxCommandEvent& event );
  void OnCheckHeatScaleClearHigh      ( wxCommandEvent& event );
  void OnCheckHeatScaleTruncate       ( wxCommandEvent& event );
  void OnCheckHeatScaleInvert         ( wxCommandEvent& event );
  
  void OnCheckDisplayVector           ( wxCommandEvent& event );
  void OnCheckDisplayTensor           ( wxCommandEvent& event );
  void OnChoiceInversion              ( wxCommandEvent& event );
  void OnChoiceRepresentation         ( wxCommandEvent& event );
  
  void OnCheckShowLabelOutline  ( wxCommandEvent& event );
  void OnChoiceUpSampleMethod   ( wxCommandEvent& event );
  
  void OnCheckContour     ( wxCommandEvent& event );
  void OnSliderContourMin ( wxScrollEvent& event );
  void OnSliderContourMax ( wxScrollEvent& event );
  void OnSliderContourMinChanging( wxScrollEvent& event );
  void OnSliderContourMaxChanging( wxScrollEvent& event );
  void OnTextContourMin   ( wxCommandEvent& event );
  void OnTextContourMax   ( wxCommandEvent& event );
  void OnCheckUseImageColorMap  ( wxCommandEvent& event );
  void OnColorContour     ( wxColourPickerEvent& event );
  void OnCheckContourExtractAll ( wxCommandEvent& event );
  void OnSliderContourSmooth    ( wxScrollEvent& event );
  void OnSliderContourSmoothChanging( wxScrollEvent& event );
  void OnTextContourSmooth      ( wxCommandEvent& event );
  void OnButtonSaveContour      ( wxCommandEvent& event );
  
  void OnCheckHideInfo    ( wxCommandEvent& event );

  void DoUpdateUI();
  void ShowWidgets( std::vector<wxWindow*>& list, bool bShow );
  void PopulateColorTable( COLOR_TABLE* ct );
  void SearchInColorTable( const wxString& search_strg );
  void UpdateColorIndicator();

  void UpdateLayerList( Layer* layer );

  virtual void DoListenToMessage( std::string const iMsg, void* iData, void* sender );

  void UpdateTextValue( wxTextCtrl* textctrl, double dvalue );

  wxCheckListBox* m_listBoxLayers;
  wxSlider*       m_sliderOpacity;
  wxTextCtrl*     m_textOpacity;
  wxCheckBox*     m_checkClearBackground;
  wxListBox*      m_listColorTable;
  wxChoice*       m_choiceColorMap;
  wxChoice*       m_choiceLUT;
  wxTextCtrl*     m_textDrawValue;
  wxTextCtrl*     m_textWindow;
  wxTextCtrl*     m_textLevel;
  wxSlider*       m_sliderWindow;
  wxSlider*       m_sliderLevel;
  wxTextCtrl*     m_textFileName;
  wxSlider*       m_sliderColorMapMin;
  wxSlider*       m_sliderColorMapMax;
  wxSlider*       m_sliderHeatScaleMid;
  wxSlider*       m_sliderHeatScaleOffset;
  wxCheckBox*     m_checkHeatScaleClearHigh;
  wxCheckBox*     m_checkHeatScaleTruncate;
  wxCheckBox*     m_checkHeatScaleInvert;
  wxTextCtrl*     m_textColorMapMin;
  wxTextCtrl*     m_textColorMapMax;
  wxTextCtrl*     m_textHeatScaleMid;
  wxTextCtrl*     m_textHeatScaleOffset;
  wxCheckBox*     m_checkSmooth;
  wxCheckBox*     m_checkUpsample;
  wxChoice*       m_choiceDirectionCode;
  wxColorIndicator* m_colorIndicator;
  wxTextCtrl*     m_textFrame;
  wxSlider*       m_sliderFrame;
  
  wxCheckBox*     m_checkDisplayVector;
  wxCheckBox*     m_checkDisplayTensor;
  wxChoice*       m_choiceInversion;
  wxChoice*       m_choiceRepresentation;
  wxChoice*       m_choiceMask;

  wxCheckBox*     m_checkShowLabelOutline;
  wxChoice*       m_choiceUpSampleMethod;
  
  wxCheckBox*     m_checkContour;
  wxSlider*       m_sliderContourMin;
  wxSlider*       m_sliderContourMax;
  wxTextCtrl*     m_textContourMin;
  wxTextCtrl*     m_textContourMax;
  wxCheckBox*     m_checkUseImageColorMap;
  wxCheckBox*     m_checkContourExtractAll;
  wxColourPickerCtrl* m_colorpickerContour;
  wxSlider*       m_sliderContourSmooth;
  wxTextCtrl*     m_textContourSmooth;
  wxButton*       m_btnSaveContour;
  
  wxCheckBox*     m_checkHideInfo;

  std::vector<wxWindow*> m_widgetlistGrayScale;
  std::vector<wxWindow*> m_widgetlistHeatScale;
  std::vector<wxWindow*> m_widgetlistGenericColorMap;
  std::vector<wxWindow*> m_widgetlistLUT;
  std::vector<wxWindow*> m_widgetlistDirectionCode;
  std::vector<wxWindow*> m_widgetlistFrame;
  std::vector<wxWindow*> m_widgetlistVector;
  std::vector<wxWindow*> m_widgetlistContour;
  std::vector<wxWindow*> m_widgetlistNormalDisplay;
  std::vector<wxWindow*> m_widgetlistEditable;
  
  std::vector<wxWindow*> m_widgetlistResize;

  LUTDataHolder* m_luts;

  COLOR_TABLE*  m_curCTAB;
  bool          m_bColorTableNeedReset;
      
  bool          m_bUINeedUpdate;

  LayerMRI*     m_layerCopied;

  DECLARE_EVENT_TABLE()
};

#endif

