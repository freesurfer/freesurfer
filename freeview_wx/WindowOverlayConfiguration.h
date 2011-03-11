/**
 * @file  WindowOverlayConfiguration.h
 * @brief Overlay configuration window.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: krish $
 *    $Date: 2011/03/11 23:27:43 $
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

#ifndef WindowOverlayConfiguration_h
#define WindowOverlayConfiguration_h

#include <wx/frame.h>
#include "Listener.h"

class wxHistogramWidget;
class wxChoice;
class wxButton;
class wxSlider;
class wxTextCtrl;
class wxRadioButton;
class wxCheckBox;
class LayerSurface;
class SurfaceOverlayProperties;

class WindowOverlayConfiguration : public wxFrame, public Listener
{
public:
  WindowOverlayConfiguration( wxWindow* parent );
	virtual ~WindowOverlayConfiguration() {}
	 
	void ShowWindow( LayerSurface* layer );

protected:
  void OnInternalIdle();
  
  void OnClose( wxCloseEvent& event );
  void OnChoiceOverlay    ( wxCommandEvent& event );
  
  void OnButtonClose      ( wxCommandEvent& event ) { Close(); }
  void OnButtonApply      ( wxCommandEvent& event );
  void OnSliderOpacity    ( wxScrollEvent& event );
  void OnTextOpacity      ( wxCommandEvent& event );
  
  /*
  void OnRadioColorScaleGreenRed    ( wxCommandEvent& event );
  void OnRadioColorScaleHeat        ( wxCommandEvent& event );
  void OnRadioColorScaleBlueRed     ( wxCommandEvent& event );
  void OnRadioColorScaleColorWheel  ( wxCommandEvent& event );
  void OnRadioColorScaleRYGBWheel   ( wxCommandEvent& event );
  void OnRadioColorScaleTwoCondGR   ( wxCommandEvent& event );
  
  void OnTextThresholdMin           ( wxCommandEvent& event );
  void OnTextThresholdMid           ( wxCommandEvent& event );
  void OnTextThresholdMax           ( wxCommandEvent& event );
  
  void OnRadioThresholdLinear       ( wxCommandEvent& event );
 */  
  
  void OnUpdateThreshold            ( wxCommandEvent& event );
  
  void OnTextThresholdChanged       ( wxCommandEvent& event );
  void OnUpdateGraph                ( wxCommandEvent& event );
  
  void UpdateGraph();
  
  void UpdateTextValue( wxTextCtrl* ctrl, double dvalue );
  
	virtual void DoListenToMessage ( std::string const iMsg, void* iData, void* sender );
  
  void UpdateUI( bool bForce = false );
  
  void DoUpdateUI();
  
  bool UpdateOverlayProperties( SurfaceOverlayProperties* p );
  
  void UpdateThresholdChanges();
			
	wxHistogramWidget*  m_histogramWidget;
	wxChoice*           m_choiceOverlay;
	
  wxSlider*           m_sliderOpacity;
  wxTextCtrl*         m_textOpacity;
  wxRadioButton*      m_radioGreenRed;
  wxRadioButton*      m_radioHeat;
  wxRadioButton*      m_radioBlueRed;
  wxRadioButton*      m_radioColorWheel;
  wxRadioButton*      m_radioRYGBWheel;
  wxRadioButton*      m_radioTwoCondGR;
  
  wxRadioButton*      m_radioThresholdLinear;
  wxRadioButton*      m_radioThresholdLinearOpaque;
  wxRadioButton*      m_radioThresholdPiecewise;
  
  wxCheckBox*         m_checkColorInverse;
  wxCheckBox*         m_checkColorTruncate;
  
  wxTextCtrl*         m_textThresholdMin;
  wxTextCtrl*         m_textThresholdMid;
  wxTextCtrl*         m_textThresholdMax;
  
  LayerSurface*			  m_layerSurface;
  
  bool                m_bUINeedUpdate; 
  
    // any class wishing to process wxWindows events must use this macro
	DECLARE_EVENT_TABLE()
};

#endif


