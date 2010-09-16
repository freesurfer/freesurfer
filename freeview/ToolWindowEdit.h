/**
 * @file  ToolWindowEdit.h
 * @brief Preferences Dialog.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2010/09/16 17:25:05 $
 *    $Revision: 1.11 $
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
#ifndef ToolWindowEdit_h
#define ToolWindowEdit_h

#include <wx/wx.h>
#include <vector>
#include "Listener.h"

class wxTextCtrl;
class wxCheckBox;
class wxComboBox;
class wxChoice;
class wxToolBar;
class wxSpinEvent;
class wxSpinCtrl;
class wxColourPickerCtrl;
class wxColourPickerEvent;

class ToolWindowEdit : public wxFrame, public Listener
{
public:
  ToolWindowEdit( wxWindow* parent );
  virtual ~ToolWindowEdit();

  void OnClose( wxCloseEvent& event );
  void OnShow( wxShowEvent& event );

  void OnActionVoxelFreehand          ( wxCommandEvent& event );
  void OnActionVoxelFreehandUpdateUI  ( wxUpdateUIEvent& event );
  void OnActionVoxelFill              ( wxCommandEvent& event );
  void OnActionVoxelFillUpdateUI      ( wxUpdateUIEvent& event );
  void OnActionVoxelPolyline          ( wxCommandEvent& event );
  void OnActionVoxelPolylineUpdateUI  ( wxUpdateUIEvent& event );
  void OnActionVoxelLivewire          ( wxCommandEvent& event );
  void OnActionVoxelLivewireUpdateUI  ( wxUpdateUIEvent& event );
  void OnActionVoxelColorPicker       ( wxCommandEvent& event );
  void OnActionVoxelColorPickerUpdateUI  ( wxUpdateUIEvent& event );
  void OnActionVoxelContour           ( wxCommandEvent& event );
  void OnActionVoxelContourUpdateUI   ( wxUpdateUIEvent& event );

  void OnActionROIFreehand            ( wxCommandEvent& event );
  void OnActionROIFreehandUpdateUI    ( wxUpdateUIEvent& event );
  void OnActionROIFill                ( wxCommandEvent& event );
  void OnActionROIFillUpdateUI        ( wxUpdateUIEvent& event );
  void OnActionROIPolyline            ( wxCommandEvent& event );
  void OnActionROIPolylineUpdateUI    ( wxUpdateUIEvent& event );
  void OnActionROILivewire            ( wxCommandEvent& event );
  void OnActionROILivewireUpdateUI    ( wxUpdateUIEvent& event );

  void OnSpinBrushSize          ( wxSpinEvent& event );
  void OnSpinBrushTolerance     ( wxSpinEvent& event );
  void OnChoiceBrushTemplate    ( wxCommandEvent& event );
  void OnCheckDrawRange         ( wxCommandEvent& event );
  void OnCheckExcludeRange      ( wxCommandEvent& event );
  void OnEditDrawRangeLow       ( wxCommandEvent& event );
  void OnEditDrawRangeHigh      ( wxCommandEvent& event );
  void OnEditExcludeRangeLow    ( wxCommandEvent& event );
  void OnEditExcludeRangeHigh   ( wxCommandEvent& event );
  void OnEditSmoothSD           ( wxCommandEvent& event );
  
  void OnEditContourValue       ( wxCommandEvent& event );
  void OnColorContour           ( wxColourPickerEvent& event );

  void OnCheckDrawConnectedOnly ( wxCommandEvent& event );
  void OnCheckSmooth            ( wxCommandEvent& event );

  void UpdateTools();

  void ResetPosition();
  
  void ShowWidgets( std::vector<wxWindow*>& list, bool bShow );

protected:
  void DoListenToMessage ( std::string const iMsg, void* iData, void* sender );
  void DoUpdateTools();
  void UpdateTextValue( wxTextCtrl* ctrl, double dvalue );

  void OnInternalIdle();

  wxToolBar*    m_toolbarVoxelEdit;
  wxToolBar*    m_toolbarROIEdit;
  wxSpinCtrl*   m_spinBrushSize;
  wxSpinCtrl*   m_spinBrushTolerance;
  wxChoice*     m_choiceTemplate;
  wxCheckBox*   m_checkDrawRange;
  wxCheckBox*   m_checkExcludeRange;
  wxTextCtrl*   m_editDrawRangeLow;
  wxTextCtrl*   m_editDrawRangeHigh;
  wxTextCtrl*   m_editExcludeRangeLow;
  wxTextCtrl*   m_editExcludeRangeHigh;
  wxCheckBox*   m_checkDrawConnectedOnly;
  wxTextCtrl*   m_editSmoothSD;
  wxCheckBox*   m_checkSmooth;
  wxTextCtrl*   m_editContourValue;
  wxColourPickerCtrl*  m_colorPickerContour;

  bool m_bToUpdateTools;
  
  std::vector<wxWindow*>  m_widgetsBrushSize;
  std::vector<wxWindow*>  m_widgetsReference;
  std::vector<wxWindow*>  m_widgetsTolerance;
  std::vector<wxWindow*>  m_widgetsConstrain;
  std::vector<wxWindow*>  m_widgetsSmooth;
  std::vector<wxWindow*>  m_widgetsContour;

  DECLARE_EVENT_TABLE()
};

#endif

