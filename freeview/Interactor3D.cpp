/**
 * @file  Interactor3D.cpp
 * @brief Interactor3D to manage mouse and key input on render view.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2008/03/27 20:38:59 $
 *    $Revision: 1.2 $
 *
 * Copyright (C) 2002-2007,
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

#include "Interactor3D.h"
#include "RenderView3D.h"
#include "MainWindow.h"
#include "LayerCollection.h"
#include "LayerCollectionManager.h"
#include "LayerPropertiesMRI.h"
#include "LayerMRI.h"
#include <vtkRenderer.h>

Interactor3D::Interactor3D() : Interactor()
{
	m_bWindowLevel = false;
}

Interactor3D::~Interactor3D()
{
}

bool Interactor3D::ProcessMouseDownEvent( wxMouseEvent& event, RenderView* renderview )
{
	m_nMousePosX = event.GetX();
	m_nMousePosY = event.GetY();
	
	if ( event.MiddleDown() && event.ControlDown() )
	{
		m_bWindowLevel = true;	
	}	
	else
	{				
		return Interactor::ProcessMouseDownEvent( event, renderview );	// pass down the event
	}	
		
	return false;	// do not pass down the event
}

bool Interactor3D::ProcessMouseUpEvent( wxMouseEvent& event, RenderView* renderview )
{	
	m_nMousePosX = event.GetX();
	m_nMousePosY = event.GetY();
	m_bWindowLevel = false;
	
	return Interactor::ProcessMouseUpEvent( event, renderview );
}

bool Interactor3D::ProcessMouseMoveEvent( wxMouseEvent& event, RenderView* renderview )
{
	LayerCollectionManager* lcm = MainWindow::GetMainWindowPointer()->GetLayerCollectionManager();
	if ( !lcm->HasAnyLayer() )
	{
		return Interactor::ProcessMouseMoveEvent( event, renderview );
	}
	
	int posX = event.GetX();
	int posY = event.GetY();
			
	if ( m_bWindowLevel )
	{
		LayerMRI* layer = ( LayerMRI* )(lcm->GetLayerCollection( "MRI" )->GetActiveLayer());
		if ( layer && layer->IsVisible() )
		{
			double scaleX = 0.002;
			double scaleY = 0.002;
			double w = ( posX - m_nMousePosX ) * scaleX;
			double l = ( posY - m_nMousePosY ) * scaleY;
			double scaleOverall = layer->GetProperties()->GetMaxValue() - 
						layer->GetProperties()->GetMinValue();
			w *= scaleOverall;
			l *= scaleOverall;
			w += layer->GetProperties()->GetWindow();
			l += layer->GetProperties()->GetLevel();
			if ( w < 0 )
				w = 0;
			layer->GetProperties()->SetWindowLevel( w, l );
		}
		m_nMousePosX = posX;
		m_nMousePosY = posY;	
	}
	else
	{	
		return Interactor::ProcessMouseMoveEvent( event, renderview );
	}
	
	return false;
}

bool Interactor3D::ProcessKeyDownEvent( wxKeyEvent& event, RenderView* renderview )
{	
	return Interactor::ProcessKeyDownEvent( event, renderview );
}
