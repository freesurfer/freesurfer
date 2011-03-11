/**
 * @file  WindowConnectivityConfiguration.h
 * @brief Connectivity display configuration window.
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

#ifndef WindowConnectivityConfiguration_h
#define WindowConnectivityConfiguration_h

#include <wx/frame.h>
#include "Listener.h"

class wxChoice;
class wxButton;
class wxSlider;
class wxTextCtrl;
class wxRadioButton;
class wxCheckBox;
class ConnectivityData;

class WindowConnectivityConfiguration : public wxFrame, public Listener
{
public:
  WindowConnectivityConfiguration( wxWindow* parent );
	virtual ~WindowConnectivityConfiguration() {}
	 
	void ShowWindow();

protected:
  void OnInternalIdle();
  
  void OnClose            ( wxCloseEvent& event ); 
  
  void OnChoiceDisplay    ( wxCommandEvent& event );  
  void OnChoiceRadius     ( wxCommandEvent& event );
  void OnTextRadiusMin    ( wxCommandEvent& event );
  void OnTextRadiusMax    ( wxCommandEvent& event );
  void OnSliderRadiusMin  ( wxScrollEvent& event );
  void OnSliderRadiusMax  ( wxScrollEvent& event );
  void OnTextThresholdMin ( wxCommandEvent& event );
  void OnTextThresholdMax ( wxCommandEvent& event );
  void OnSliderThresholdMin ( wxScrollEvent& event );
  void OnSliderThresholdMax ( wxScrollEvent& event );
  void OnCheckAddOn       ( wxCommandEvent& event );
  void OnButtonExport     ( wxCommandEvent& event );
  
  void UpdateTextValue( wxTextCtrl* ctrl, double dvalue );
  
	virtual void DoListenToMessage ( std::string const iMsg, void* iData, void* sender );
  
  void UpdateUI( bool bForce = false );
  
  void DoUpdateUI();
 
  wxChoice*           m_choiceDisplay;	
  wxChoice*           m_choiceRadius;
  wxSlider*           m_sliderRadiusMin;
  wxSlider*           m_sliderRadiusMax;
  wxTextCtrl*         m_textRadiusMin;
  wxTextCtrl*         m_textRadiusMax;
  wxSlider*           m_sliderThresholdMin;
  wxSlider*           m_sliderThresholdMax;
  wxTextCtrl*         m_textThresholdMin;
  wxTextCtrl*         m_textThresholdMax;
  wxCheckBox*         m_checkAddOn;
  wxButton*           m_btnExport;

  std::vector<wxWindow*>  m_widgetsRadiusMax;
  
  ConnectivityData*	  m_conn;
  
  bool                m_bUINeedUpdate; 
  
    // any class wishing to process wxWindows events must use this macro
	DECLARE_EVENT_TABLE()
};

#endif


