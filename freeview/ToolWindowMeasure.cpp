/**
 * @file  ToolWindowMeasure.h
 * @brief Preferences dialog.
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



#include "ToolWindowMeasure.h"
#include <wx/wx.h>
#include <wx/config.h>
#include <wx/xrc/xmlres.h>
#include <wx/filedlg.h>
#include <wx/filename.h>
#include <wx/spinctrl.h>
#include <wx/listbox.h>
#include <wx/clipbrd.h>
#include <wx/ffile.h>
#include "MainWindow.h"
#include "RenderView2D.h"
#include "RenderView3D.h"
#include "BrushProperty.h"
#include "Interactor2DMeasure.h"
#include "LayerCollection.h"
#include "LayerMRI.h"
#include "LayerPropertiesMRI.h"
#include "Region2D.h"

BEGIN_EVENT_TABLE( ToolWindowMeasure, wxFrame )
  EVT_MENU      ( XRCID( "ID_ACTION_MEASURE_LINE" ),    ToolWindowMeasure::OnActionMeasureLine )
  EVT_UPDATE_UI ( XRCID( "ID_ACTION_MEASURE_LINE" ),    ToolWindowMeasure::OnActionMeasureLineUpdateUI )
  EVT_MENU      ( XRCID( "ID_ACTION_MEASURE_RECT" ),    ToolWindowMeasure::OnActionMeasureRectangle )
  EVT_UPDATE_UI ( XRCID( "ID_ACTION_MEASURE_RECT" ),    ToolWindowMeasure::OnActionMeasureRectangleUpdateUI )

  EVT_BUTTON    ( XRCID( "ID_BUTTON_COPY" ),            ToolWindowMeasure::OnButtonCopy )
  EVT_BUTTON    ( XRCID( "ID_BUTTON_EXPORT" ),          ToolWindowMeasure::OnButtonExport )
  
  EVT_SHOW      ( ToolWindowMeasure::OnShow )

END_EVENT_TABLE()


ToolWindowMeasure::ToolWindowMeasure( wxWindow* parent ) : Listener( "ToolWindowMeasure" )
{
  wxXmlResource::Get()->LoadFrame( this, parent, wxT("ID_TOOLWINDOW_MEASURE") );
  m_toolbar   = XRCCTRL( *this, "ID_TOOLBAR_MEASURE", wxToolBar );
  m_textStats = XRCCTRL( *this, "ID_TEXT_STATS",      wxTextCtrl );
  m_btnCopy   = XRCCTRL( *this, "ID_BUTTON_COPY",     wxButton );
  m_btnExport = XRCCTRL( *this, "ID_BUTTON_EXPORT",   wxButton );
  m_region = NULL;
  m_bToUpdateStats = false;
}

ToolWindowMeasure::~ToolWindowMeasure()
{}

void ToolWindowMeasure::OnShow( wxShowEvent& event )
{
  if ( event.GetShow() )
  {
    wxConfigBase* config = wxConfigBase::Get();
    if ( config )
    {
      int x = config->Read( _T("/ToolWindowMeasure/PosX"), 50L );
      int y = config->Read( _T("/ToolWindowMeasure/PosY"), 50L );
      Move( x, y );
    }
    
 //   SetClientSize( GetClientSize().GetWidth(), m_toolbar->GetSize().GetHeight() );
  }
  else
  {
    wxConfigBase* config = wxConfigBase::Get();
    if ( config )
    {
      int x, y;
      GetPosition( &x, &y );
      config->Write( _T("/ToolWindowMeasure/PosX"), (long) x );
      config->Write( _T("/ToolWindowMeasure/PosY"), (long) y );
    }
  }
}

void ToolWindowMeasure::ResetPosition()
{
  if ( IsShown() )
  {
    wxConfigBase* config = wxConfigBase::Get();
    if ( config )
    {
      int x = config->Read( _T("/ToolWindowMeasure/PosX"), 50L );
      int y = config->Read( _T("/ToolWindowMeasure/PosY"), 50L );
      Move( x, y );
    }
  }
}

void ToolWindowMeasure::SetRegion( Region2D* reg )
{
  m_region = reg;
  if ( m_region )
    m_region->AddListener( this );  
  UpdateStats();
}

void ToolWindowMeasure::UpdateStats( )
{
  m_bToUpdateStats = true;
}

void ToolWindowMeasure::DoUpdateStats()
{
  wxString strg;
  if ( m_region )
  {
    wxArrayString strgs = m_region->GetLongStats();
    for ( size_t i = 0; i < strgs.size(); i++ )
      strg += strgs[i] + "\n";   
  }
  m_textStats->ChangeValue( strg ); 
  m_btnCopy->Enable( !strg.IsEmpty() );
  m_btnExport->Enable( !strg.IsEmpty() );
  
  m_bToUpdateStats = false;
}

void ToolWindowMeasure::OnInternalIdle()
{
  wxFrame::OnInternalIdle();
  
  if ( m_bToUpdateStats )
    DoUpdateStats();
}

void ToolWindowMeasure::DoListenToMessage ( std::string const iMsg, void* iData, void* sender )
{
  if ( iMsg == "RegionStatsUpdated" )
  {
    if ( m_region == iData || m_region == sender )
      UpdateStats();
  }
}

void ToolWindowMeasure::OnActionMeasureLine( wxCommandEvent& event )
{
  MainWindow::GetMainWindowPointer()->SetAction( Interactor2DMeasure::MM_Line );
}

void ToolWindowMeasure::OnActionMeasureLineUpdateUI( wxUpdateUIEvent& event)
{
  RenderView2D* view = ( RenderView2D* )MainWindow::GetMainWindowPointer()->GetRenderView( 0 );
  event.Check( view->GetInteractionMode() == RenderView2D::IM_Measure
      && view->GetAction() == Interactor2DMeasure::MM_Line );
  event.Enable( view->GetInteractionMode() == RenderView2D::IM_Measure
                && !MainWindow::GetMainWindowPointer()->GetLayerCollection( "MRI" )->IsEmpty() );
}


void ToolWindowMeasure::OnActionMeasureRectangle( wxCommandEvent& event )
{
  MainWindow::GetMainWindowPointer()->SetAction( Interactor2DMeasure::MM_Rectangle );
}

void ToolWindowMeasure::OnActionMeasureRectangleUpdateUI( wxUpdateUIEvent& event)
{
  RenderView2D* view = ( RenderView2D* )MainWindow::GetMainWindowPointer()->GetRenderView( 0 );
  event.Check( view->GetInteractionMode() == RenderView2D::IM_Measure
      && view->GetAction() == Interactor2DMeasure::MM_Rectangle );
  event.Enable( view->GetInteractionMode() == RenderView2D::IM_Measure
      && !MainWindow::GetMainWindowPointer()->GetLayerCollection( "MRI" )->IsEmpty() );
}

void ToolWindowMeasure::OnButtonCopy( wxCommandEvent& event )
{
  wxString output = m_textStats->GetValue();
  if (wxTheClipboard->Open())
  {
    wxTheClipboard->SetData( new wxTextDataObject( output ) );
    wxTheClipboard->Close();
  }
}

void ToolWindowMeasure::OnButtonExport( wxCommandEvent& event )
{
  wxFileDialog dlg( this, _("Export to file"), _(""), _(""),
                    _("All files (*.*)|*.*"),
                    wxFD_SAVE | wxFD_OVERWRITE_PROMPT );
  if ( dlg.ShowModal() == wxID_OK )
  {
    wxString fn = dlg.GetPath();
    if ( !m_textStats->SaveFile( fn ) )
    {
      wxMessageDialog msg_dlg( this, wxString("Can not write to file ") + fn, 
                           _("Error"), wxOK );
      msg_dlg.ShowModal();
    } 
  }
}

