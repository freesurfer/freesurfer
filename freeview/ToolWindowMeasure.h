/**
 * @file  ToolWindowMeasure.h
 * @brief Preferences Dialog.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2009/11/03 22:51:29 $
 *    $Revision: 1.1 $
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

class wxToolBar;

class ToolWindowMeasure : public wxFrame
{
public:
  ToolWindowMeasure( wxWindow* parent );
  virtual ~ToolWindowMeasure();

  void OnShow( wxShowEvent& event );

  void OnActionMeasureLine              ( wxCommandEvent& event );
  void OnActionMeasureLineUpdateUI      ( wxUpdateUIEvent& event );
  void OnActionMeasureRectangle         ( wxCommandEvent& event );
  void OnActionMeasureRectangleUpdateUI ( wxUpdateUIEvent& event );

  void ResetPosition();

protected:
  wxToolBar*  m_toolbar;
  
  DECLARE_EVENT_TABLE()
};

#endif

