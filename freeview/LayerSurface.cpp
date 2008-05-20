/**
 * @file  LayerSurface.cpp
 * @brief Layer data object for MRI surface.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2008/05/20 16:28:32 $
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
#include "vtkRGBATransferFunction.h"
#include "vtkLookupTable.h"
#include "vtkProperty.h"
#include "vtkImageReslice.h"
#include "vtkFreesurferLookupTable.h"
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
	}
	mProperties = new LayerPropertiesSurface();
	mProperties->AddListener( this );
	m_mainActor = vtkActor::New();
}

LayerSurface::~LayerSurface()
{	
	for ( int i = 0; i < 3; i++ )
	{
		m_sliceActor2D[i]->Delete();
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
	
	event.SetInt( 100 );
	wxPostEvent( wnd, event );
	
//	InitializeSurface();
//	InitializeActors();
	
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
/*	m_dWorldOrigin[0] = RASBounds[0];
	m_dWorldOrigin[1] = RASBounds[2];
	m_dWorldOrigin[2] = RASBounds[4]; 
	source->GetPixelSize( m_dWorldVoxelSize );
	
	m_dWorldSize[0] = ( ( int )( (RASBounds[1] - RASBounds[0]) / m_dWorldVoxelSize[0] ) ) * m_dWorldVoxelSize[0];
	m_dWorldSize[1] = ( ( int )( (RASBounds[3] - RASBounds[2]) / m_dWorldVoxelSize[1] ) ) * m_dWorldVoxelSize[1];
	m_dWorldSize[2] = ( ( int )( (RASBounds[5] - RASBounds[4]) / m_dWorldVoxelSize[2] ) ) * m_dWorldVoxelSize[2];
	
	m_volumeRAS = source->GetImageOutput();	*/
}
	
		
void LayerSurface::InitializeActors()
{
/*	for ( int i = 0; i < 3; i++ )
	{
  // The reslice object just takes a slice out of the volume.
		//
		mReslice[i] = vtkSmartPointer<vtkImageReslice>::New();
		mReslice[i]->SetInput( m_volumeRAS );
//		mReslice[i]->SetOutputSpacing( sizeX, sizeY, sizeZ );
		mReslice[i]->BorderOff();

  // This sets us to extract slices.
		mReslice[i]->SetOutputDimensionality( 2 );

  // This will change depending what orienation we're in.
		mReslice[i]->SetResliceAxesDirectionCosines( 1, 0, 0,
				0, 1, 0,
				0, 0, 1 );

  // This will change to select a different slice.
		mReslice[i]->SetResliceAxesOrigin( 0, 0, 0 );
				
		//
  // Image to colors using color table.
		//
		mColorMap[i] = vtkSmartPointer<vtkImageMapToColors>::New();
		mColorMap[i]->SetLookupTable( mProperties->GetGrayScaleTable() );
		mColorMap[i]->SetInput( mReslice[i]->GetOutput() );
		mColorMap[i]->SetOutputFormatToRGBA();
		mColorMap[i]->PassAlphaToOutputOn();

		//
  // Prop in scene with plane mesh and texture.
		//
		m_sliceActor2D[i]->SetInput( mColorMap[i]->GetOutput() );
		m_sliceActor3D[i]->SetInput( mColorMap[i]->GetOutput() );

  // Set ourselves up.
		this->OnSlicePositionChanged( i );		
	}
		
	this->UpdateResliceInterpolation();
	this->UpdateTextureSmoothing();
	this->UpdateOpacity();
	this->UpdateColorMap();*/
}

void LayerSurface::UpdateOpacity()
{
	for ( int i = 0; i < 3; i++ )
	{
	//	m_sliceActor2D[i]->SetOpacity( mProperties->GetOpacity() );
	//	m_sliceActor3D[i]->SetOpacity( mProperties->GetOpacity() );
	}	
}

void LayerSurface::Append2DProps( vtkRenderer* renderer, int nPlane )
{
	wxASSERT ( nPlane >= 0 && nPlane <= 2 );
	
	renderer->AddViewProp( m_sliceActor2D[nPlane] );
}

void LayerSurface::Append3DProps( vtkRenderer* renderer )
{
	renderer->AddViewProp( m_mainActor ); 
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
	
	/*	
	vtkSmartPointer<vtkMatrix4x4> matrix = 
			vtkSmartPointer<vtkMatrix4x4>::New();
	matrix->Identity();
	switch ( nPlane ) 
	{
		case 0:		
			m_sliceActor2D[0]->PokeMatrix( matrix );
			m_sliceActor2D[0]->SetPosition( m_dSlicePosition[0], 0, 0 );
			m_sliceActor2D[0]->RotateX( 90 );
			m_sliceActor2D[0]->RotateY( -90 );
			m_sliceActor3D[0]->PokeMatrix( matrix );
			m_sliceActor3D[0]->SetPosition( m_dSlicePosition[0], 0, 0 );
			m_sliceActor3D[0]->RotateX( 90 );
			m_sliceActor3D[0]->RotateY( -90 );
    
    // Putting negatives in the reslice axes cosines will flip the
    // image on that axis.
			mReslice[0]->SetResliceAxesDirectionCosines( 0, -1, 0,
					0, 0, 1,
					1, 0, 0 );
			mReslice[0]->SetResliceAxesOrigin( m_dSlicePosition[0], 0, 0  );
			
			break;
		case 1:
			m_sliceActor2D[1]->PokeMatrix( matrix );
			m_sliceActor2D[1]->SetPosition( 0, m_dSlicePosition[1], 0 );
			m_sliceActor2D[1]->RotateX( 90 );
		//	m_sliceActor2D[1]->RotateY( 180 );
			m_sliceActor3D[1]->PokeMatrix( matrix );
			m_sliceActor3D[1]->SetPosition( 0, m_dSlicePosition[1], 0 );
			m_sliceActor3D[1]->RotateX( 90 );
		//	m_sliceActor3D[1]->RotateY( 180 );
    
    // Putting negatives in the reslice axes cosines will flip the
      // image on that axis.
			mReslice[1]->SetResliceAxesDirectionCosines( 1, 0, 0,
					0, 0, 1,
					0, 1, 0 );
			mReslice[1]->SetResliceAxesOrigin( 0, m_dSlicePosition[1], 0 );
			break;
		case 2:
			m_sliceActor2D[2]->SetPosition( 0, 0, m_dSlicePosition[2] );
		//	m_sliceActor2D[2]->RotateY( 180 );
			m_sliceActor3D[2]->SetPosition( 0, 0, m_dSlicePosition[2] );
		//	m_sliceActor3D[2]->RotateY( 180 );
			
			mReslice[2]->SetResliceAxesDirectionCosines( 1, 0, 0,
					0, 1, 0,
					0, 0, 1 );
			mReslice[2]->SetResliceAxesOrigin( 0, 0, m_dSlicePosition[2]  );
			break;
	}
	*/
}

LayerPropertiesSurface* LayerSurface::GetProperties()
{
	return mProperties;
}

void LayerSurface::DoListenToMessage( std::string const iMessage, void* const iData )
{
/*	if ( iMessage == "ColorMapChanged" )
	{
		this->UpdateColorMap();
		this->SendBroadcast( "LayerActorUpdated", this );
	}
	else if ( iMessage == "ResliceInterpolationChanged" )
	{
		this->UpdateResliceInterpolation();
		this->SendBroadcast( "LayerActorUpdated", this );
	}
	else if ( iMessage == "OpacityChanged" )
	{
		this->UpdateOpacity();
		this->SendBroadcast( "LayerActorUpdated", this );
	}
	else if ( iMessage == "TextureSmoothingChanged" )
	{
		this->UpdateTextureSmoothing();
		this->SendBroadcast( "LayerActorUpdated", this );
	}
	else if ( iMessage == "WindowLevelChanged" )
	{
		this->SendBroadcast( iMessage, this );
}*/
}

void LayerSurface::SetVisible( bool bVisible )
{
	for ( int i = 0; i < 3; i++ )
	{
		m_sliceActor2D[i]->SetVisibility( bVisible ? 1 : 0 );
	}
	m_mainActor->SetVisibility( bVisible ? 1 : 0 );
	this->SendBroadcast( "LayerActorUpdated", this );
}

bool LayerSurface::IsVisible()
{
	return m_sliceActor2D[0]->GetVisibility() > 0;
}

