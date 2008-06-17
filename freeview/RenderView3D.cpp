/**
 * @file  RenderView3D.cpp
 * @brief View for rendering.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2008/06/17 23:08:18 $
 *    $Revision: 1.4 $
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
 
#include "RenderView3D.h"
#include "MainWindow.h"
#include "LayerCollection.h"
#include "LayerMRI.h"
#include <vtkRenderer.h>
#include "vtkConeSource.h"
#include "vtkPolyDataMapper.h"
#include "vtkActor.h"
#include "vtkCamera.h"
#include "vtkActor2D.h"
#include "vtkCellPicker.h"
#include "vtkPointPicker.h"
#include "vtkProp3DCollection.h"
#include "Interactor3DNavigate.h"

IMPLEMENT_DYNAMIC_CLASS(RenderView3D, RenderView)

BEGIN_EVENT_TABLE(RenderView3D, RenderView)

END_EVENT_TABLE()

RenderView3D::RenderView3D() : RenderView()
{
	InitializeRenderView3D();
}

RenderView3D::RenderView3D( wxWindow* parent, int id ) : RenderView( parent, id )
{
	InitializeRenderView3D();
}

void RenderView3D::InitializeRenderView3D()
{
	this->SetDesiredUpdateRate( 5000 );
//	this->SetStillUpdateRate( 0.5 );
		
	if ( m_interactor )
		delete m_interactor;
	
	m_interactor = new Interactor3DNavigate();
	
	m_bToUpdateRASPosition = false;
	vtkCellPicker* picker = vtkCellPicker::New();
//	vtkPointPicker* picker = vtkPointPicker::New();
	this->SetPicker( picker );
	picker->Delete();
}

RenderView3D* RenderView3D::New()
{
  // we don't make use of the objectfactory, because we're not registered
  return new RenderView3D;
}

RenderView3D::~RenderView3D()
{

}

void RenderView3D::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}

void RenderView3D::RefreshAllActors()
{
	LayerCollectionManager* lcm = MainWindow::GetMainWindowPointer()->GetLayerCollectionManager();
	
	m_renderer->RemoveAllViewProps();
	lcm->Append3DProps( m_renderer );
	
	// add focus frame
	m_renderer->AddViewProp( m_actorFocusFrame );
	
	m_renderer->ResetCameraClippingRange();
	
	Render();
}

void RenderView3D::UpdateViewByWorldCoordinate()
{
	vtkCamera* cam = m_renderer->GetActiveCamera();
	double wcenter[3];
	for ( int i = 0; i < 3; i++ )
	{
		wcenter[i] = m_dWorldOrigin[i] + m_dWorldSize[i] / 2;
	}
	cam->SetFocalPoint( wcenter );
	cam->SetPosition( wcenter[0]+m_dWorldSize[0]*2.5, wcenter[1] + m_dWorldSize[1]*.25, wcenter[2]+m_dWorldSize[2]*.25 );
	cam->SetViewUp( 0, 0, 1 );
	m_renderer->ResetCameraClippingRange();
}

void RenderView3D::UpdateMouseRASPosition( int posX, int posY )
{
	m_bToUpdateRASPosition = true;
	m_nPickCoord[0] = posX;
	m_nPickCoord[1] = posY;
}

void RenderView3D::CancelUpdateMouseRASPosition()
{
	m_bToUpdateRASPosition = false;
}

void RenderView3D::DoUpdateMouseRASPosition( int posX, int posY )
{	
	LayerCollection* lc_mri = MainWindow::GetMainWindowPointer()->GetLayerCollection( "MRI" );	
	LayerCollection* lc_roi = MainWindow::GetMainWindowPointer()->GetLayerCollection( "ROI" );	
	LayerCollection* lc_surface = MainWindow::GetMainWindowPointer()->GetLayerCollection( "Surface" );	
	
//	MousePositionToRAS( posX, posY, pos );
//	vtkPointPicker* picker = vtkPointPicker::SafeDownCast( this->GetPicker() );
	vtkCellPicker* picker = vtkCellPicker::SafeDownCast( this->GetPicker() );
	if ( picker )
	{
		double pos[3];
		picker->Pick( posX, GetClientSize().GetHeight() - posY, 0, GetRenderer() );
		picker->GetPickPosition( pos );

		vtkProp* prop = picker->GetViewProp();
	//	cout << pos[0] << " " << pos[1] << " " << pos[2] << ",   " << prop << endl;
		if ( prop && ( lc_mri->HasProp( prop ) || lc_roi->HasProp( prop ) || lc_surface->HasProp( prop ) ) )
		{		
			lc_mri->SetCurrentRASPosition( pos );
		}
	}
}

void RenderView3D::OnInternalIdle()
{
	RenderView::OnInternalIdle();
	
	if ( m_bToUpdateRASPosition )
	{
		DoUpdateMouseRASPosition( m_nPickCoord[0], m_nPickCoord[1] );
		m_bToUpdateRASPosition = false;
	}
		
}
