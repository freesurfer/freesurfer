/**
 * @file  WindowHistogram.h
 * @brief Main window.
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

#ifndef WindowHistogram_h
#define WindowHistogram_h

#include <wx/frame.h>
#include "Listener.h"

class wxHistogramWidget;
class wxChoice;
class Layer;

class WindowHistogram : public wxFrame, public Listener
{
public:
  WindowHistogram( wxWindow* parent );
	virtual ~WindowHistogram() {}
	
	void ShowWindow( Layer* default_layer );
	void UpdateHistogram();

protected:
  void OnClose( wxCloseEvent& event );
  void OnChoiceLayer( wxCommandEvent& event );
  void OnButtonClose( wxCommandEvent& event ) { Close(); }
    
	void DoListenToMessage ( std::string const iMsg, void* iData, void* sender );
	void UpdateUI();
			
	wxHistogramWidget*	m_histogramWidget;
	wxChoice*			m_choiceLayer;
	
	Layer*			m_activeLayer;
    // any class wishing to process wxWindows events must use this macro
	DECLARE_EVENT_TABLE()
};

#endif


