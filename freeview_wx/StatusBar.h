/**
 * @file  StatusBar.h
 * @brief Main window StatusBar.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: krish $
 *    $Date: 2011/03/11 23:27:42 $
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

#ifndef StatusBar_h
#define StatusBar_h

#include "wx/statusbr.h"

class wxGauge;
class wxTimer;
class wxTimerEvent;

class StatusBar : public wxStatusBar
{
public:
  StatusBar( wxWindow* parent );
  virtual ~StatusBar();

  void OnSize( wxSizeEvent& event );
  void OnTimer( wxTimerEvent& event );

  void ActivateProgressBar();

  wxGauge* m_gaugeBar;

protected:
  wxTimer* m_timer;

  // any class wishing to process wxWindows events must use this macro
  DECLARE_EVENT_TABLE()

};

#endif


