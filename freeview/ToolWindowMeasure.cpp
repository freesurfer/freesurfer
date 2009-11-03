/**
 * @file  ToolWindowMeasure.h
 * @brief Preferences dialog.
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



#include "ToolWindowMeasure.h"
#include <wx/wx.h>
#include <wx/config.h>
#include <wx/xrc/xmlres.h>
#include <wx/filedlg.h>
#include <wx/filename.h>
#include <wx/spinctrl.h>
#include "MainWindow.h"
#include "RenderView2D.h"
#include "RenderView3D.h"
#include "BrushProperty.h"
#include "Interactor2DMeasure.h"
#include "LayerCollection.h"
#include "LayerMRI.h"
#include "LayerPropertiesMRI.h"

BEGIN_EVENT_TABLE( ToolWindowMeasure, wxFrame )
EVT_MENU      ( XRCID( "ID_ACTION_MEASURE_LINE" ),    ToolWindowMeasure::OnActionMeasureLine )
EVT_UPDATE_UI ( XRCID( "ID_ACTION_MEASURE_LINE" ),    ToolWindowMeasure::OnActionMeasureLineUpdateUI )
EVT_MENU      ( XRCID( "ID_ACTION_MEASURE_RECT" ),    ToolWindowMeasure::OnActionMeasureRectangle )
EVT_UPDATE_UI ( XRCID( "ID_ACTION_MEASURE_RECT" ),    ToolWindowMeasure::OnActionMeasureRectangleUpdateUI )

EVT_SHOW      ( ToolWindowMeasure::OnShow )

END_EVENT_TABLE()


ToolWindowMeasure::ToolWindowMeasure( wxWindow* parent ) 
{
  wxXmlResource::Get()->LoadFrame( this, parent, wxT("ID_TOOLWINDOW_MEASURE") );
  m_toolbar = XRCCTRL( *this, "ID_TOOLBAR_MEASURE", wxToolBar );
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
