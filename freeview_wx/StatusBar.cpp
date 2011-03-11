/**
 * @file  StatusBar.cpp
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

#include "StatusBar.h"
#include <wx/gauge.h>
#include <wx/timer.h>

#define PROGRESS_TIMER_ID 0

BEGIN_EVENT_TABLE( StatusBar, wxStatusBar )
  EVT_SIZE  ( StatusBar::OnSize )
  EVT_TIMER ( PROGRESS_TIMER_ID, StatusBar::OnTimer )
END_EVENT_TABLE()


StatusBar::StatusBar( wxWindow* parent ) : wxStatusBar( parent )
{
  static const int widths[2] =
    {
      -1, 200
    };
  static const int styles[2] =
    {
      wxSB_NORMAL, wxSB_NORMAL
    };
  SetFieldsCount( 2 );
  SetStatusWidths( 2, widths );
  SetStatusStyles( 2, styles );

  m_gaugeBar = new wxGauge( this, wxID_ANY, 100 );
  m_gaugeBar->Hide();
  m_gaugeBar->SetValue(0);

  m_timer = new wxTimer( this, PROGRESS_TIMER_ID );
}

StatusBar::~StatusBar()
{}

void StatusBar::OnSize( wxSizeEvent& event )
{
  wxRect rect;
  GetFieldRect( 1, rect );
  m_gaugeBar->SetSize( rect.x, rect.y, rect.width, rect.height );
  event.Skip();
}

void StatusBar::OnTimer( wxTimerEvent& event )
{
  if ( m_gaugeBar->IsShown() )
  {
    int val = m_gaugeBar->GetValue() + 1;
    if (val > m_gaugeBar->GetRange())
      val = m_gaugeBar->GetRange();
    m_gaugeBar->SetValue( val );
  }
  else
    m_timer->Stop();
}

void StatusBar::ActivateProgressBar()
{
  m_gaugeBar->Show();
  m_gaugeBar->SetValue( 0 );

  m_timer->Stop();
  m_timer->Start( 250 );
}
