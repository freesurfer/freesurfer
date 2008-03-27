/**
 * @file  Annotation2D.cpp
 * @brief View for rendering.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2008/03/27 20:38:58 $
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
 
#include <wx/xrc/xmlres.h> 
#include "Annotation2D.h"
#include "vtkRenderer.h"
#include "vtkTextActor.h"
#include "vtkTextProperty.h"
#include "LayerCollection.h"
#include "MainWindow.h"
#include "vtkPropCollection.h"
#include "LayerMRI.h"
#include "LayerPropertiesMRI.h"

#define NUMBER_OF_ANNOTATIONS	5


Annotation2D::Annotation2D()
{
	for ( int i = 0; i < NUMBER_OF_ANNOTATIONS; i++ )
	{
		m_actorCoordinates[i] = vtkSmartPointer<vtkTextActor>::New();
		m_actorCoordinates[i]->GetPositionCoordinate()->SetCoordinateSystemToNormalizedViewport();
		m_actorCoordinates[i]->GetTextProperty()->ShadowOff();
		m_actorCoordinates[i]->GetTextProperty()->SetFontSize(11);
		m_actorCoordinates[i]->GetTextProperty()->ItalicOff();
	//	m_actorCoordinates[i]->GetTextProperty()->SetFontFamilyToTimes();
	}
	m_actorCoordinates[0]->SetPosition( 0.01, 0.5 );
	m_actorCoordinates[0]->GetTextProperty()->SetJustificationToLeft();
	m_actorCoordinates[0]->GetTextProperty()->SetVerticalJustificationToCentered();
	m_actorCoordinates[1]->SetPosition( 0.5, 0.99 );
	m_actorCoordinates[1]->GetTextProperty()->SetJustificationToCentered();
	m_actorCoordinates[1]->GetTextProperty()->SetVerticalJustificationToTop();
	m_actorCoordinates[2]->SetPosition( 0.99, 0.5 );
	m_actorCoordinates[2]->GetTextProperty()->SetJustificationToRight();
	m_actorCoordinates[2]->GetTextProperty()->SetVerticalJustificationToCentered();
	m_actorCoordinates[3]->SetPosition( 0.5, 0.01 );
	m_actorCoordinates[3]->GetTextProperty()->SetJustificationToCentered();
	m_actorCoordinates[3]->GetTextProperty()->SetVerticalJustificationToBottom();
	
	// indicate slice location
	m_actorCoordinates[4]->SetPosition( 0.99, 0.01 );
	m_actorCoordinates[4]->GetTextProperty()->SetJustificationToRight();
	m_actorCoordinates[4]->GetTextProperty()->SetVerticalJustificationToBottom();
	
}

Annotation2D::~Annotation2D()
{
}

void Annotation2D::Update( vtkRenderer* renderer, int nPlane )
{
	double slicePos[3] = { 0, 0, 0 };	
	LayerCollection* lc = MainWindow::GetMainWindowPointer()->GetLayerCollection( "MRI" );
	lc->GetSlicePosition( slicePos );
	bool bHasLayer = ( lc->GetNumberOfLayers() > 0 );
	for ( int i = 0; i < NUMBER_OF_ANNOTATIONS; i++ )
		m_actorCoordinates[i]->SetVisibility( bHasLayer ? 1 : 0 );
			
	if ( !bHasLayer )
	{
		return;
	}
	
	LayerMRI* mri = ( LayerMRI* )lc->GetLayer( 0 );
	double centPos[3];
	double* worigin = mri->GetWorldOrigin();
	double* wsize = mri->GetWorldSize();
	for ( int i = 0; i < 3; i++ )
		centPos[i] = worigin[i] + wsize[i]/2;
	
	double pos[4] = { 0, 1, 1, 0 }, tpos = 0;
	
	switch ( nPlane )
	{
		case 0:
			renderer->NormalizedViewportToView( pos[0], pos[1], tpos );
			renderer->ViewToWorld( pos[0], pos[1], tpos );
			pos[0] = pos[1]; pos[1] = tpos;
			renderer->NormalizedViewportToView( pos[2], pos[3], tpos );
			renderer->ViewToWorld( pos[2], pos[3], tpos );
			pos[2] = pos[3]; pos[3] = tpos;
			
			tpos = slicePos[0];
			mri->RemapPositionToRealRAS( tpos, pos[0], pos[1], slicePos[0], pos[0], pos[1] );
			if ( pos[0] >= 0 )
				m_actorCoordinates[0]->SetInput( wxString::Format( "A %.2f", pos[0]).c_str() );
			else
				m_actorCoordinates[0]->SetInput( wxString::Format( "P %.2f", -pos[0]).c_str() );
			
			if ( pos[1] >= 0 )
				m_actorCoordinates[1]->SetInput( wxString::Format( "S %.2f", pos[1]).c_str() );
			else
				m_actorCoordinates[1]->SetInput( wxString::Format( "I  %.2f", -pos[1]).c_str() );
			
			mri->RemapPositionToRealRAS( tpos, pos[2], pos[3], slicePos[0], pos[2], pos[3] );
			if ( pos[2] >= 0 )
				m_actorCoordinates[2]->SetInput( wxString::Format( "A %.2f", pos[2]).c_str() );
			else
				m_actorCoordinates[2]->SetInput( wxString::Format( "P %.2f", -pos[2]).c_str() );
			
			if ( pos[3] >= 0 )
				m_actorCoordinates[3]->SetInput( wxString::Format( "S %.2f", pos[3]).c_str() );
			else
				m_actorCoordinates[3]->SetInput( wxString::Format( "I  %.2f", -pos[3]).c_str() );
						
			mri->RemapPositionToRealRAS( tpos, centPos[1], centPos[2], slicePos[0], centPos[1], centPos[2] );
			if ( slicePos[0] >= 0 )
				m_actorCoordinates[4]->SetInput( wxString::Format( "R %.2f", slicePos[0]).c_str() );
			else
				m_actorCoordinates[4]->SetInput( wxString::Format( "L %.2f", -slicePos[0]).c_str() );
			
			break;
		case 1:
			renderer->NormalizedViewportToView( pos[0], pos[1], tpos );
			renderer->ViewToWorld( pos[0], pos[1], tpos );
			pos[1] = tpos;
			renderer->NormalizedViewportToView( pos[2], pos[3], tpos );
			renderer->ViewToWorld( pos[2], pos[3], tpos );
			pos[3] = tpos;	
			
			tpos = slicePos[1];
			mri->RemapPositionToRealRAS( pos[0], tpos, pos[1], pos[0], slicePos[1], pos[1] );
			if ( pos[0] >= 0 )
				m_actorCoordinates[0]->SetInput( wxString::Format( "R %.2f", pos[0]).c_str() );
			else 
				m_actorCoordinates[0]->SetInput( wxString::Format( "L %.2f", -pos[0]).c_str() );
			
			if ( pos[1] >= 0 )
				m_actorCoordinates[1]->SetInput( wxString::Format( "S %.2f", pos[1]).c_str() );
			else
				m_actorCoordinates[1]->SetInput( wxString::Format( "I  %.2f", -pos[1]).c_str() );
			
			mri->RemapPositionToRealRAS( pos[2], tpos, pos[3], pos[2], slicePos[1], pos[3] );
			
			if ( pos[2] >= 0 )
				m_actorCoordinates[2]->SetInput( wxString::Format( "R %.2f", pos[2]).c_str() );
			else
				m_actorCoordinates[2]->SetInput( wxString::Format( "L %.2f", -pos[2]).c_str() );
			
			if ( pos[3] >= 0 )
				m_actorCoordinates[3]->SetInput( wxString::Format( "S %.2f", pos[3]).c_str() );
			else
				m_actorCoordinates[3]->SetInput( wxString::Format( "I  %.2f", -pos[3]).c_str() );
			
			mri->RemapPositionToRealRAS( centPos[0], tpos, centPos[2], centPos[0], slicePos[1], centPos[2] );
			if ( slicePos[1] >= 0 )
				m_actorCoordinates[4]->SetInput( wxString::Format( "A %.2f", slicePos[1]).c_str() );
			else
				m_actorCoordinates[4]->SetInput( wxString::Format( "P %.2f", -slicePos[1]).c_str() );
			
			break;
		case 2:
			renderer->NormalizedViewportToView( pos[0], pos[1], tpos );
			renderer->ViewToWorld( pos[0], pos[1], tpos );
			tpos = 0;
			renderer->NormalizedViewportToView( pos[2], pos[3], tpos );
			renderer->ViewToWorld( pos[2], pos[3], tpos );			
			
			tpos = slicePos[2];
			mri->RemapPositionToRealRAS( pos[0], pos[1], tpos, pos[0], pos[1], slicePos[2] );

			if ( pos[0] >= 0 )
				m_actorCoordinates[0]->SetInput( wxString::Format( "R %.2f", pos[0]).c_str() );
			else 
				m_actorCoordinates[0]->SetInput( wxString::Format( "L %.2f", -pos[0]).c_str() );
			
			if ( pos[1] >= 0 )
				m_actorCoordinates[1]->SetInput( wxString::Format( "A %.2f", pos[1]).c_str() );
			else
				m_actorCoordinates[1]->SetInput( wxString::Format( "P %.2f", -pos[1]).c_str() );
			
			mri->RemapPositionToRealRAS( pos[2], pos[3], tpos, pos[2], pos[3], slicePos[2] );
			
			if ( pos[2] >= 0 )
				m_actorCoordinates[2]->SetInput( wxString::Format( "R %.2f", pos[2]).c_str() );
			else
				m_actorCoordinates[2]->SetInput( wxString::Format( "L %.2f", -pos[2]).c_str() );
			
			if ( pos[3] >= 0 )
				m_actorCoordinates[3]->SetInput( wxString::Format( "A %.2f", pos[3]).c_str() );
			else
				m_actorCoordinates[3]->SetInput( wxString::Format( "P %.2f", -pos[3]).c_str() );
			
			mri->RemapPositionToRealRAS( centPos[0], centPos[1], tpos, centPos[0], centPos[1], slicePos[2] );
			if ( slicePos[2] >= 0 )
				m_actorCoordinates[4]->SetInput( wxString::Format( "S %.2f", slicePos[2]).c_str() );
			else
				m_actorCoordinates[4]->SetInput( wxString::Format( "I  %.2f", -slicePos[2]).c_str() );
			
			break;
	}
}

void Annotation2D::AppendAnnotations( vtkRenderer* renderer )
{
	for ( int i = 0; i < NUMBER_OF_ANNOTATIONS; i++ )
		renderer->AddViewProp( m_actorCoordinates[i] );
}
