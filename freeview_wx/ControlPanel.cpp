/**
 * @file  ControlPanel.h
 * @brief Main control panel.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: krish $
 *    $Date: 2011/03/11 23:27:36 $
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

#include "wx/wx.h"

#include "ControlPanel.h"
#include "PanelVolume.h"
#include "PanelROI.h"
#include "PanelSurface.h"
#include "PanelWayPoints.h"
#include "PanelSceneSetting.h"
#include <wx/xrc/xmlres.h>
#include <wx/aui/auibook.h>
#include <wx/laywin.h>

BEGIN_EVENT_TABLE( ControlPanel, wxPanel )
/* EVT_IDLE(ControlPanel::OnIdle)
 EVT_CHECKBOX  (XRCID(wxT("ID_CHECKBOX_SLICE_X")),   ControlPanel::OnCheckBoxSlice)
 EVT_CHECKBOX  (XRCID(wxT("ID_CHECKBOX_SLICE_Y")),   ControlPanel::OnCheckBoxSlice)
 EVT_CHECKBOX  (XRCID(wxT("ID_CHECKBOX_SLICE_Z")),   ControlPanel::OnCheckBoxSlice)
 EVT_TEXT   (XRCID(wxT("ID_TEXTCTRL_SLICE_X")),   ControlPanel::OnEditSlice)
 EVT_TEXT   (XRCID(wxT("ID_TEXTCTRL_SLICE_Y")),   ControlPanel::OnEditSlice)
 EVT_TEXT   (XRCID(wxT("ID_TEXTCTRL_SLICE_Z")),   ControlPanel::OnEditSlice)
 EVT_COMMAND_SCROLL (XRCID(wxT("ID_SLIDER_SLICE_X")),   ControlPanel::OnSliderSlice)
 EVT_COMMAND_SCROLL (XRCID(wxT("ID_SLIDER_SLICE_Y")),   ControlPanel::OnSliderSlice)
 EVT_COMMAND_SCROLL (XRCID(wxT("ID_SLIDER_SLICE_Z")),   ControlPanel::OnSliderSlice)
 EVT_CHECKBOX  (XRCID(wxT("ID_CHECKBOX_ATTACH_SLICE")),  ControlPanel::OnCheckBoxAttachSlice)
 EVT_COMMAND_SCROLL (XRCID(wxT("ID_SLIDER_SLICE_TRANSPARENCY")), ControlPanel::OnSliderTransparency)*/
END_EVENT_TABLE()


ControlPanel::ControlPanel( wxWindow* parent ) : wxPanel( parent, wxID_ANY)
{
  wxBoxSizer* sizer = new wxBoxSizer( wxVERTICAL );

  m_notebook = new wxAuiNotebook( this );
  m_notebook->SetWindowStyle( wxAUI_NB_TOP | wxAUI_NB_TAB_SPLIT | wxAUI_NB_TAB_MOVE |
                              wxAUI_NB_SCROLL_BUTTONS | wxAUI_NB_WINDOWLIST_BUTTON );
  m_panelVolume = new PanelVolume( m_notebook );
  m_notebook->AddPage( m_panelVolume, _("Volumes") );

  m_panelSurface = new PanelSurface( m_notebook );
  m_notebook->AddPage( m_panelSurface, _("Surfaces") );

  m_panelROI = new PanelROI( m_notebook );
  m_notebook->AddPage( m_panelROI, _("ROIs") );
  
  m_panelWayPoints = new PanelWayPoints( m_notebook );
  m_notebook->AddPage( m_panelWayPoints, _("Point Sets") );

// m_panelSceneSetting = new PanelSceneSetting( m_notebook );
// m_notebook->AddPage( m_panelSceneSetting, _T("Scene Settings") );

  sizer->Add( m_notebook, 1, wxEXPAND | wxALL );

  this->SetSizer( sizer );
  this->Layout();
}

ControlPanel::~ControlPanel()
{}

void ControlPanel::RaisePage( const wxString& title )
{
  for ( unsigned int i = 0; i < m_notebook->GetPageCount(); i++ )
  {
    if ( m_notebook->GetPageText( i ) == title )
    {
      m_notebook->SetSelection( i );
      return;
    }
  }
}

void ControlPanel::UpdateUI( bool bForce )
{
  m_panelVolume->UpdateUI( bForce );
  m_panelSurface->UpdateUI( bForce );
  m_panelROI->UpdateUI( bForce );
  m_panelWayPoints->UpdateUI( bForce );
}

wxString ControlPanel::GetCurrentLayerCollectionName()
{
  return m_notebook->GetPageText( m_notebook->GetSelection() );
}

