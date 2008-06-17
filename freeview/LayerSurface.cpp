/**
 * @file  LayerSurface.cpp
 * @brief Layer data object for MRI surface.
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
 
#include <wx/wx.h>
#include "LayerSurface.h"
#include "vtkRenderer.h"
#include "vtkImageActor.h"
#include "vtkPolyDataMapper.h"
#include "vtkSmartPointer.h"
#include "vtkMatrix4x4.h"
#include "vtkImageMapToColors.h"
#include "vtkImageActor.h"
#include "vtkActor.h"
#include "vtkLODActor.h"
#include "vtkRGBATransferFunction.h"
#include "vtkLookupTable.h"
#include "vtkProperty.h"
#include "vtkImageReslice.h"
#include "vtkFreesurferLookupTable.h"
#include "vtkPlane.h"
#include "vtkCutter.h"
#include "vtkDecimatePro.h"
#include "LayerPropertiesSurface.h"
#include "MyUtils.h"
#include "FSSurface.h"

LayerSurface::LayerSurface() : Layer(),
		m_surfaceSource( NULL),
		m_bResampleToRAS( true )
{
	m_strTypeNames.push_back( "Surface" );
	
	for ( int i = 0; i < 3; i++ )
	{
	//	m_nSliceNumber[i] = 0;
		m_sliceActor2D[i] = vtkActor::New();
		m_sliceActor3D[i] = vtkActor::New();
	}
	mProperties = new LayerPropertiesSurface();
	mProperties->AddListener( this );
	m_mainActor = vtkLODActor::New();
	mLowResFilter = vtkSmartPointer<vtkDecimatePro>::New();
	mLowResFilter->SetTargetReduction( 0.9 );
//	mMediumResFilter = vtkSmartPointer<vtkDecimatePro>::New();
//	mMediumResFilter->SetTargetReduction( 0.9 );
}

LayerSurface::~LayerSurface()
{	
	for ( int i = 0; i < 3; i++ )
	{
		m_sliceActor2D[i]->Delete();
		m_sliceActor3D[i]->Delete();
	}
	m_mainActor->Delete();
	
	delete mProperties;
	
	if ( m_surfaceSource )
		delete m_surfaceSource;
}
 
bool LayerSurface::LoadSurfaceFromFile( wxWindow* wnd, wxCommandEvent& event )
{
	if ( m_surfaceSource )
		delete m_surfaceSource;
	
	m_surfaceSource = new FSSurface();
//	m_surfaceSource->SetResampleToRAS( m_bResampleToRAS );
	if ( !m_surfaceSource->MRISRead( m_sFilename.c_str(), wnd, event ) )
		return false;
	
	InitializeSurface();
	InitializeActors();
	
	event.SetInt( 100 );
	wxPostEvent( wnd, event );
	
//	mProperties->SetSurfaceSource( m_surfaceSource );
	
	return true;	
}

void LayerSurface::InitializeSurface()
{
	if ( m_surfaceSource == NULL )
		return;

	FSSurface* source = m_surfaceSource;
		
	float RASBounds[6];
	source->GetBounds( RASBounds );	
	m_dWorldOrigin[0] = RASBounds[0];
	m_dWorldOrigin[1] = RASBounds[2];
	m_dWorldOrigin[2] = RASBounds[4]; 
	
	m_dWorldSize[0] = RASBounds[1] - RASBounds[0];
	m_dWorldSize[1] = RASBounds[3] - RASBounds[2];
	m_dWorldSize[2] = RASBounds[5] - RASBounds[4];
}
	
		
void LayerSurface::InitializeActors()
{
	if ( m_surfaceSource == NULL )
		return;

	vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	mapper->SetInput( m_surfaceSource->GetPolyData() );
	m_mainActor->SetMapper( mapper );
	mapper->Update();
	mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	mLowResFilter->SetInput( m_surfaceSource->GetPolyData() );
	mapper->SetInput( mLowResFilter->GetOutput() );
	m_mainActor->AddLODMapper( mapper );
	mapper->Update();
/*	mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	mMediumResFilter->SetInput( m_surfaceSource->GetPolyData() );
	mapper->SetInput( mMediumResFilter->GetOutput() );
	m_mainActor->AddLODMapper( mapper );
	mapper->Update();*/
	
	for ( int i = 0; i < 3; i++ )
	{
  // The reslice object just takes a slice out of the volume.
		//
		mReslicePlane[i] = vtkSmartPointer<vtkPlane>::New();
		mReslicePlane[i]->SetOrigin( 0, 0, 0 );
		mReslicePlane[i]->SetNormal( (i==0), (i==1), (i==2) );

		vtkSmartPointer<vtkCutter> cutter = 
				vtkSmartPointer<vtkCutter>::New();
		cutter->SetInput( m_surfaceSource->GetPolyData() );
		cutter->SetCutFunction( mReslicePlane[i] );
		    
		//
    // Mappers for the lines.
		//
		vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
		mapper->SetInput( cutter->GetOutput() );
		vtkSmartPointer<vtkPolyDataMapper> mapper2 = vtkSmartPointer<vtkPolyDataMapper>::New();
		mapper2->SetInput( cutter->GetOutput() );
		//
    // Actors in the scene, drawing the mapped lines.
		//
		m_sliceActor2D[i]->SetMapper( mapper );
//		m_sliceActor2D[i]->SetBackfaceProperty( m_sliceActor2D[i]->MakeProperty() );
//		m_sliceActor2D[i]->GetBackfaceProperty()->BackfaceCullingOff();
		m_sliceActor2D[i]->SetProperty( m_sliceActor2D[i]->MakeProperty() );
		m_sliceActor2D[i]->GetProperty()->SetInterpolationToFlat();
		m_sliceActor2D[i]->GetProperty()->SetLineWidth( 1 );
		
		m_sliceActor3D[i]->SetMapper( mapper2 );
//		m_sliceActor3D[i]->SetBackfaceProperty( m_sliceActor3D[i]->MakeProperty() );
//		m_sliceActor3D[i]->GetBackfaceProperty()->BackfaceCullingOff();
		m_sliceActor3D[i]->SetProperty( m_sliceActor3D[i]->MakeProperty() );
		m_sliceActor3D[i]->GetProperty()->SetLineWidth( 2 );
		m_sliceActor3D[i]->GetProperty()->SetInterpolationToFlat();		

  // Set ourselves up.
		this->OnSlicePositionChanged( i );	
	}
		
	this->UpdateOpacity();
	this->UpdateColorMap();
}

void LayerSurface::UpdateOpacity()
{
	for ( int i = 0; i < 3; i++ )
	{
	//	m_sliceActor2D[i]->GetProperty()->SetOpacity( mProperties->GetOpacity() );
	//	m_sliceActor3D[i]->SetOpacity( mProperties->GetOpacity() );
	}	
	m_mainActor->GetProperty()->SetOpacity( mProperties->GetOpacity() );
}

void LayerSurface::UpdateColorMap () 
{
	assert( mProperties );

	for ( int i = 0; i < 3; i++ )
	{
		m_sliceActor2D[i]->GetProperty()->SetColor( mProperties->GetEdgeColor() );
		m_sliceActor3D[i]->GetProperty()->SetColor( mProperties->GetEdgeColor() );
	}
	
	m_mainActor->GetProperty()->SetColor( mProperties->GetColor() );
}

void LayerSurface::Append2DProps( vtkRenderer* renderer, int nPlane )
{
	wxASSERT ( nPlane >= 0 && nPlane <= 2 );
	
	renderer->AddViewProp( m_sliceActor2D[nPlane] );
}

void LayerSurface::Append3DProps( vtkRenderer* renderer )
{
	renderer->AddViewProp( m_mainActor ); 
	
	for ( int i = 0; i < 3; i++ )
		renderer->AddViewProp( m_sliceActor3D[i] );
}

/*
void LayerSurface::SetSliceNumber( int* sliceNumber )
{
	if ( sliceNumber[0] != m_nSliceNumber[0] || sliceNumber[1] != m_nSliceNumber[1] ||
			sliceNumber[2] != m_nSliceNumber[2] )
	{
		m_nSliceNumber[0] = sliceNumber[0];
		m_nSliceNumber[1] = sliceNumber[1];
		m_nSliceNumber[2] = sliceNumber[2];
		
		
	}
}
*/

void LayerSurface::SetSlicePositionToWorldCenter()
{
	if ( m_surfaceSource == NULL )
		return;
	
	// Get some values from the MRI.
	double pos[3];
	for ( int i = 0; i < 3; i++ )
		pos[i] = ((int)( m_dWorldSize[i]/2/m_dWorldVoxelSize[i] ) + 0.0 ) * m_dWorldVoxelSize[i] + m_dWorldOrigin[i];
	
	SetSlicePosition( pos );
}

void LayerSurface::OnSlicePositionChanged( int nPlane ) 
{  
	if ( m_surfaceSource == NULL )
		return;
	
	assert( mProperties );
		
	vtkSmartPointer<vtkMatrix4x4> matrix = 
			vtkSmartPointer<vtkMatrix4x4>::New();
	matrix->Identity();
	switch ( nPlane ) 
	{
		case 0:		
			mReslicePlane[0]->SetOrigin( m_dSlicePosition[0], 0, 0  );	
			m_sliceActor2D[0]->SetPosition(	0.1, 0, 0 );	
			break;
		case 1:
			mReslicePlane[1]->SetOrigin( 0, m_dSlicePosition[1], 0 );
			m_sliceActor2D[1]->SetPosition(	0, 0.1, 0 );	
			break;
		case 2:
			mReslicePlane[2]->SetOrigin( 0, 0, m_dSlicePosition[2]  );
			m_sliceActor2D[2]->SetPosition(	0, 0, -0.1 );	
			break;
	}
}

LayerPropertiesSurface* LayerSurface::GetProperties()
{
	return mProperties;
}

void LayerSurface::DoListenToMessage( std::string const iMessage, void* const iData )
{
	if ( iMessage == "ColorMapChanged" )
	{
		this->UpdateColorMap();
		this->SendBroadcast( "LayerActorUpdated", this );
	}
	else if ( iMessage == "OpacityChanged" )
	{
		this->UpdateOpacity();
		this->SendBroadcast( "LayerActorUpdated", this );
	}
}

void LayerSurface::SetVisible( bool bVisible )
{
	for ( int i = 0; i < 3; i++ )
	{
		m_sliceActor2D[i]->SetVisibility( bVisible ? 1 : 0 );
		m_sliceActor3D[i]->SetVisibility( bVisible ? 1 : 0 );
	}
	m_mainActor->SetVisibility( bVisible ? 1 : 0 );
	this->SendBroadcast( "LayerActorUpdated", this );
}

bool LayerSurface::IsVisible()
{
	return m_sliceActor2D[0]->GetVisibility() > 0;
}

bool LayerSurface::HasProp( vtkProp* prop )
{
	for ( int i = 0; i < 3; i++ )
	{
		if ( m_sliceActor2D[i] == prop || m_sliceActor3D[i] == prop )
			return true;
	}
	return prop == m_mainActor;
}

int LayerSurface::GetVertexIndexAtRAS( double* ras, double* distance )
{
	if ( m_surfaceSource == NULL )
		return -1;
	
	return m_surfaceSource->FindVertexAtRAS( ras, distance );
}
