/**
 * @file  ToolWindowMeasure.h
 * @brief Preferences Dialog.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2010/02/19 16:33:44 $
 *    $Revision: 1.3 $
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
#ifndef ToolWindowMeasure_h
#define ToolWindowMeasure_h

#include <wx/wx.h>
#include "Listener.h"

class wxToolBar;
class wxTextCtrl;
class wxButton;
class Region2D;

class ToolWindowMeasure : public wxFrame, public Listener
{
public:
  ToolWindowMeasure( wxWindow* parent );
  virtual ~ToolWindowMeasure();

  void OnShow( wxShowEvent& event );

  void OnActionMeasureLine              ( wxCommandEvent& event );
  void OnActionMeasureLineUpdateUI      ( wxUpdateUIEvent& event );
  void OnActionMeasureRectangle         ( wxCommandEvent& event );
  void OnActionMeasureRectangleUpdateUI ( wxUpdateUIEvent& event );
  
  void OnButtonCopy         ( wxCommandEvent& event );
  void OnButtonExport       ( wxCommandEvent& event );

  void ResetPosition();
  
  void UpdateStats();
  
  void SetRegion( Region2D* reg );

protected:
  void DoListenToMessage ( std::string const iMsg, void* iData, void* sender );
  
  void OnInternalIdle();
 
  void DoUpdateStats();
  
  wxToolBar*  m_toolbar;
  wxTextCtrl* m_textStats;
  wxButton*   m_btnExport;
  wxButton*   m_btnCopy;
  
  Region2D*   m_region;
  bool        m_bToUpdateStats;
  
  DECLARE_EVENT_TABLE()
};

#endif

