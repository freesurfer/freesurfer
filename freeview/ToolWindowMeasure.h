/**
 * @file  ToolWindowMeasure.h
 * @brief Preferences Dialog.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2010/06/22 20:48:31 $
 *    $Revision: 1.6 $
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
class wxSpinCtrl;
class wxSpinEvent;
class wxButton;
class Region2D;
class SurfaceRegion;

class ToolWindowMeasure : public wxFrame, public Listener
{
public:
  ToolWindowMeasure( wxWindow* parent );
  virtual ~ToolWindowMeasure();

  void OnShow( wxShowEvent& event );

  void OnActionMeasureLine                  ( wxCommandEvent& event );
  void OnActionMeasureLineUpdateUI          ( wxUpdateUIEvent& event );
  void OnActionMeasureRectangle             ( wxCommandEvent& event );
  void OnActionMeasureRectangleUpdateUI     ( wxUpdateUIEvent& event );
  void OnActionMeasurePolyline              ( wxCommandEvent& event );
  void OnActionMeasurePolylineUpdateUI      ( wxUpdateUIEvent& event );
  void OnActionMeasureSpline                ( wxCommandEvent& event );
  void OnActionMeasureSplineUpdateUI        ( wxUpdateUIEvent& event );
  void OnActionMeasureSurfaceRegion         ( wxCommandEvent& event );
  void OnActionMeasureSurfaceRegionUpdateUI ( wxUpdateUIEvent& event );
  
  void OnButtonCopy         ( wxCommandEvent& event );
  void OnButtonExport       ( wxCommandEvent& event );
  void OnButtonSave         ( wxCommandEvent& event );
  void OnButtonSaveAll      ( wxCommandEvent& event );
  void OnButtonLoad         ( wxCommandEvent& event );
  void OnSpinId             ( wxSpinEvent& evnet );
  
  void UpdateWidgets();
  
  void SetRegion( Region2D* reg );
  
  void SetSurfaceRegion( SurfaceRegion* reg );

protected:
  void DoListenToMessage ( std::string const iMsg, void* iData, void* sender );
  
  void OnInternalIdle();
 
  void DoUpdateWidgets();
  
  wxToolBar*      m_toolbar;
  wxTextCtrl*     m_textStats;
  wxButton*       m_btnExport;
  wxButton*       m_btnCopy;
  wxButton*       m_btnSave;
  wxButton*       m_btnSaveAll;
  wxButton*       m_btnLoad;
  wxSpinCtrl*     m_spinId;
  
  std::vector<wxWindow*>  m_widgets2D;    // widgets for 2D measurements
  std::vector<wxWindow*>  m_widgets3D;    // widgets for 3D measurements
  
  Region2D*       m_region;
  SurfaceRegion*  m_surfaceRegion;
  bool            m_bToUpdateWidgets;
  
  DECLARE_EVENT_TABLE()
};

#endif

