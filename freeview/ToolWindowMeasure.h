/**
 * @file  ToolWindowMeasure.h
 * @brief Preferences Dialog.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 22:00:37 $
 *    $Revision: 1.11 $
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
#ifndef ToolWindowMeasure_h
#define ToolWindowMeasure_h

#include <wx/wx.h>
#include "Listener.h"

class wxToolBar;
class wxTextCtrl;
class wxSpinCtrl;
class wxSpinEvent;
class wxColourPickerCtrl;
class wxColourPickerEvent;
class wxButton;
class Region2D;
class SurfaceRegion;

class ToolWindowMeasure : public wxFrame, public Listener
{
public:
  ToolWindowMeasure( wxWindow* parent );
  virtual ~ToolWindowMeasure();

  void OnShow( wxShowEvent& event );
  void OnClose( wxCloseEvent& event);

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
  void OnActionMeasureLabel                ( wxCommandEvent& event );
  void OnActionMeasureLabelUpdateUI       ( wxUpdateUIEvent& event );
  
  void OnButtonCopy         ( wxCommandEvent& event );
  void OnButtonExport       ( wxCommandEvent& event );
  void OnButtonSave         ( wxCommandEvent& event );
  void OnButtonSaveAll      ( wxCommandEvent& event );
  void OnButtonLoad         ( wxCommandEvent& event );
  void OnButtonUpdate       ( wxCommandEvent& event );
  void OnSpinId             ( wxSpinEvent& evnet );
  void OnSpinGroup          ( wxSpinEvent& evnet );
  void OnColorGroup         ( wxColourPickerEvent& event );
  
  void UpdateWidgets();
  
  void SetRegion( Region2D* reg );
  
  void SetSurfaceRegion( SurfaceRegion* reg );

protected:
  void DoListenToMessage ( std::string const iMsg, void* iData, void* sender );
  
  void OnInternalIdle();
 
  wxString GetLabelStats();
  
  void DoUpdateWidgets();
  
  wxToolBar*      m_toolbar;
  wxTextCtrl*     m_textStats;
  wxButton*       m_btnExport;
  wxButton*       m_btnCopy;
  wxButton*       m_btnSave;
  wxButton*       m_btnSaveAll;
  wxButton*       m_btnLoad;
  wxButton*       m_btnUpdate;
  wxSpinCtrl*     m_spinId;
  wxSpinCtrl*     m_spinGroup;
  wxColourPickerCtrl* m_colorPickerGroup;
  
  std::vector<wxWindow*>  m_widgets2D;    // widgets for 2D measurements
  std::vector<wxWindow*>  m_widgets3D;    // widgets for 3D measurements
  
  Region2D*       m_region;
  SurfaceRegion*  m_surfaceRegion;
  bool            m_bToUpdateWidgets;
  
  DECLARE_EVENT_TABLE()
};

#endif

