/**
 * @file  RenderView3D.cpp
 * @brief View for rendering.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2008/03/27 18:12:15 $
 *    $Revision: 1.1 $
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
#include "Interactor3DNavigate.h"

IMPLEMENT_DYNAMIC_CLASS(RenderView3D, RenderView)

BEGIN_EVENT_TABLE(RenderView3D, RenderView)

END_EVENT_TABLE()

RenderView3D::RenderView3D() : RenderView()
{
	if ( m_interactor )
		delete m_interactor;
	
	m_interactor = new Interactor3DNavigate();
}

RenderView3D::RenderView3D( wxWindow* parent, int id ) : RenderView( parent, id )
{
	if ( m_interactor )
		delete m_interactor;
	
	m_interactor = new Interactor3DNavigate();
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
