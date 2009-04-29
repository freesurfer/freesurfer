/**
 * @file  ToolWindowEdit.h
 * @brief Preferences Dialog.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2009/04/29 22:53:55 $
 *    $Revision: 1.6.2.2 $
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

class wxTextCtrl;
class wxCheckBox;
class wxComboBox;
class wxChoice;
class wxToolBar;
class wxSpinEvent;
class wxSpinCtrl;

class ToolWindowEdit : public wxFrame
{
public:
  ToolWindowEdit( wxWindow* parent );
  virtual ~ToolWindowEdit();

  void OnShow( wxShowEvent& event );

  void OnActionVoxelFreehand( wxCommandEvent& event );
  void OnActionVoxelFreehandUpdateUI( wxUpdateUIEvent& event );
  void OnActionVoxelFill( wxCommandEvent& event );
  void OnActionVoxelFillUpdateUI( wxUpdateUIEvent& event );
  void OnActionVoxelPolyline( wxCommandEvent& event );
  void OnActionVoxelPolylineUpdateUI( wxUpdateUIEvent& event );
  void OnActionVoxelLivewire( wxCommandEvent& event );
  void OnActionVoxelLivewireUpdateUI( wxUpdateUIEvent& event );

  void OnActionROIFreehand( wxCommandEvent& event );
  void OnActionROIFreehandUpdateUI( wxUpdateUIEvent& event );
  void OnActionROIFill( wxCommandEvent& event );
  void OnActionROIFillUpdateUI( wxUpdateUIEvent& event );
  void OnActionROIPolyline( wxCommandEvent& event );
  void OnActionROIPolylineUpdateUI( wxUpdateUIEvent& event );
  void OnActionROILivewire( wxCommandEvent& event );
  void OnActionROILivewireUpdateUI( wxUpdateUIEvent& event );

  void OnSpinBrushSize( wxSpinEvent& event );
  void OnSpinBrushTolerance( wxSpinEvent& event );
  void OnChoiceBrushTemplate( wxCommandEvent& event );
  void OnCheckDrawRange( wxCommandEvent& event );
  void OnCheckExcludeRange( wxCommandEvent& event );
  void OnEditDrawRangeLow( wxCommandEvent& event );
  void OnEditDrawRangeHigh( wxCommandEvent& event );
  void OnEditExcludeRangeLow( wxCommandEvent& event );
  void OnEditExcludeRangeHigh( wxCommandEvent& event );

  void OnCheckDrawConnectedOnly( wxCommandEvent& event );

  void UpdateTools();

  void ResetPosition();

protected:
  void DoUpdateTools();
  void UpdateTextValue( wxTextCtrl* ctrl, double dvalue );

  void OnInternalIdle();

  wxToolBar*  m_toolbarVoxelEdit;
  wxToolBar*  m_toolbarROIEdit;
  wxSpinCtrl*  m_spinBrushSize;
  wxSpinCtrl*  m_spinBrushTolerance;
  wxChoice*  m_choiceTemplate;
  wxCheckBox*  m_checkDrawRange;
  wxCheckBox*  m_checkExcludeRange;
  wxTextCtrl*  m_editDrawRangeLow;
  wxTextCtrl*  m_editDrawRangeHigh;
  wxTextCtrl*  m_editExcludeRangeLow;
  wxTextCtrl*  m_editExcludeRangeHigh;
  wxCheckBox*  m_checkDrawConnectedOnly;

  bool m_bToUpdateTools;

  DECLARE_EVENT_TABLE()
};

#endif

