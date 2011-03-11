/**
 * @file  Interactor2DWayPointsEdit.cpp
 * @brief Interactor for editing way points in 2D render view.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: krish $
 *    $Date: 2011/03/11 23:27:38 $
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

#include "Interactor2DWayPointsEdit.h"
#include "RenderView2D.h"
#include "Cursor2D.h"
#include "MainWindow.h"
#include "LayerCollection.h"
#include "LayerCollectionManager.h"
#include "LayerPropertiesMRI.h"
#include "LayerWayPoints.h"
#include "LayerMRI.h"
#include "CursorFactory.h"
#include <vtkRenderer.h>

Interactor2DWayPointsEdit::Interactor2DWayPointsEdit() :
    Interactor2D(),
    m_bEditing( false ),
    m_nCurrentIndex( -1 )
{}

Interactor2DWayPointsEdit::~Interactor2DWayPointsEdit()
{}

bool Interactor2DWayPointsEdit::ProcessMouseDownEvent( wxMouseEvent& event, RenderView* renderview )
{
  RenderView2D* view = ( RenderView2D* )renderview;

  if ( event.LeftDown() )
  {
    if ( event.CmdDown() && event.ShiftDown() )
      return Interactor2D::ProcessMouseDownEvent( event, renderview );

    LayerCollection* lc = MainWindow::GetMainWindowPointer()->GetLayerCollectionManager()->GetLayerCollection( "WayPoints" );
    LayerWayPoints* wp = ( LayerWayPoints* )lc->GetActiveLayer();
    if ( !wp || !wp->IsVisible() )
    {
      SendBroadcast( "WayPointsNotVisible", this );
    }
    else
    {
      wp->SaveForUndo();
      m_nMousePosX = event.GetX();
      m_nMousePosY = event.GetY();
      m_bEditing = true;

      double ras[3];
      view->MousePositionToRAS( m_nMousePosX, m_nMousePosY, ras );
      // wp->SetVoxelByRAS( ras, view->GetViewPlane(), !event.ShiftDown() );
      if ( !event.ShiftDown() )
      {
        m_nCurrentIndex = wp->FindPoint( ras );
        if ( m_nCurrentIndex < 0 )
          m_nCurrentIndex = wp->AddPoint( ras );
      }
      else
        wp->RemovePoint( ras );
    }

    return false;
  }
  else
    return Interactor2D::ProcessMouseDownEvent( event, renderview ); // pass down the event
}

bool Interactor2DWayPointsEdit::ProcessMouseUpEvent( wxMouseEvent& event, RenderView* renderview )
{
// RenderView2D* view = ( RenderView2D* )renderview;
  UpdateCursor( event, renderview );

  m_bEditing = false;
//  m_nCurrentIndex = -1;

  return Interactor2D::ProcessMouseUpEvent( event, renderview );
}

bool Interactor2DWayPointsEdit::ProcessMouseMoveEvent( wxMouseEvent& event, RenderView* renderview )
{
  RenderView2D* view = ( RenderView2D* )renderview;			
	if ( m_bEditing )
	{
		UpdateCursor( event, view );	
		int posX = event.GetX();
		int posY = event.GetY();
		m_nMousePosX = posX;
		m_nMousePosY = posY;
		if ( m_nCurrentIndex >= 0 )
		{
			double ras[3];
			view->MousePositionToRAS( m_nMousePosX, m_nMousePosY, ras );			
			LayerCollection* lc = MainWindow::GetMainWindowPointer()->GetLayerCollection( "WayPoints" );
			LayerWayPoints* wp = ( LayerWayPoints* )lc->GetActiveLayer();
			wp->UpdatePoint( ras, m_nCurrentIndex );
		}
		
		return false;
	} 
	else
	{	
		return Interactor2D::ProcessMouseMoveEvent( event, renderview );
	}
}

bool Interactor2DWayPointsEdit::ProcessKeyDownEvent( wxKeyEvent& event, RenderView* renderview )
{
  UpdateCursor( event, renderview );
  
  if ( !m_bEditing )
    return Interactor2D::ProcessKeyDownEvent( event, renderview );
  else
    return false;
}

void Interactor2DWayPointsEdit::UpdateCursor( wxEvent& event, wxWindow* wnd )
{
  if ( wnd->FindFocus() == wnd )
  {
    if ( event.IsKindOf( CLASSINFO( wxMouseEvent ) ) )
    {
      wxMouseEvent* e = ( wxMouseEvent* )&event;
      if ( e->MiddleDown() || e->RightDown() ||
           ( e->CmdDown() && e->ShiftDown() ) )
      {
        Interactor2D::UpdateCursor( event, wnd );
      }
      else
      {
        // set own cursor
        wnd->SetCursor( wxNullCursor );
      }
    }
  }
  else
    Interactor2D::UpdateCursor( event, wnd );
}
