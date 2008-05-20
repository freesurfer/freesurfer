/**
 * @file  LayerMRI.cpp
 * @brief Layer data object for MRI volume.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2008/05/20 16:28:32 $
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
#include "LayerMRI.h"
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
#include "LayerPropertiesMRI.h"
#include "MyUtils.h"
#include "FSVolume.h"

LayerMRI::LayerMRI() : LayerEditable(),
		m_volumeSource( NULL),
		m_bResampleToRAS( true )
{
	m_strTypeNames.push_back( "MRI" );
	
	for ( int i = 0; i < 3; i++ )
	{
	//	m_nSliceNumber[i] = 0;
		m_sliceActor2D[i] = vtkImageActor::New();
		m_sliceActor3D[i] = vtkImageActor::New();
	/*	m_sliceActor2D[i]->GetProperty()->SetAmbient( 1 );
		m_sliceActor2D[i]->GetProperty()->SetDiffuse( 0 );
		m_sliceActor3D[i]->GetProperty()->SetAmbient( 1 );
		m_sliceActor3D[i]->GetProperty()->SetDiffuse( 0 );*/
		m_sliceActor2D[i]->InterpolateOff();
		m_sliceActor3D[i]->InterpolateOff();
	}
	mProperties = new LayerPropertiesMRI();
	mProperties->AddListener( this );
}

LayerMRI::~LayerMRI()
{	
	for ( int i = 0; i < 3; i++ )
	{
		m_sliceActor2D[i]->Delete();
		m_sliceActor3D[i]->Delete();
	}
	
	delete mProperties;
	
	if ( m_volumeSource )
		delete m_volumeSource;
}
 
/*
bool LayerMRI::LoadVolumeFromFile( std::string filename )
{
	m_sFilename = filename;
	
	return LoadVolumeFromFile();
}
*/

void LayerMRI::SetResampleToRAS( bool bResample )
{
	m_bResampleToRAS = bResample;
}

bool LayerMRI::LoadVolumeFromFile( wxWindow* wnd, wxCommandEvent& event )
{
	if ( m_volumeSource )
		delete m_volumeSource;
	
	m_volumeSource = new FSVolume();
	m_volumeSource->SetResampleToRAS( m_bResampleToRAS );
	if ( !m_volumeSource->MRIRead( m_sFilename.c_str(), wnd, event ) )
		return false;
	
	event.SetInt( 100 );
	wxPostEvent( wnd, event );
	
	InitializeVolume();
	InitializeActors();
	
	mProperties->SetVolumeSource( m_volumeSource );
	
	return true;	
}

bool LayerMRI::Create( LayerMRI* mri, bool bCopyVoxelData )
{
	if ( m_volumeSource )
		delete m_volumeSource;
	
	m_volumeSource = new FSVolume();
	m_volumeSource->Create( mri->m_volumeSource, bCopyVoxelData );
	m_bResampleToRAS = mri->m_bResampleToRAS;
	m_volumeRef = mri->GetRASVolume();
	if ( m_volumeRef.GetPointer() )
	{
		SetWorldOrigin( mri->GetWorldOrigin() );
		SetWorldVoxelSize( mri->GetWorldVoxelSize() );
		SetWorldSize( mri->GetWorldSize() );
		
		m_volumeRAS = vtkSmartPointer<vtkImageData>::New();
		
		if ( bCopyVoxelData )
		{
			m_volumeRAS->DeepCopy( m_volumeRef );
		}
		else
		{			
			m_volumeRAS->SetNumberOfScalarComponents( 1 );
			m_volumeRAS->SetScalarTypeToFloat();
			m_volumeRAS->SetOrigin( m_volumeRef->GetOrigin() );
			m_volumeRAS->SetSpacing( m_volumeRef->GetSpacing() );	
			m_volumeRAS->SetDimensions( m_volumeRef->GetDimensions() );
			m_volumeRAS->AllocateScalars();
			float* ptr = ( float* )m_volumeRAS->GetScalarPointer();
			int* nDim = m_volumeRAS->GetDimensions();
		//	cout << nDim[0] << ", " << nDim[1] << ", " << nDim[2] << endl;
			memset( ptr, 0, sizeof( float ) * nDim[0] * nDim[1] * nDim[2] );
		}
		
		InitializeActors();
		
		mProperties->SetVolumeSource( m_volumeSource );
	//	mProperties->SetColorMap( LayerPropertiesMRI::LUT );
		
		m_sFilename = "";
		
		if ( bCopyVoxelData )
			SetModified();
	}
	
	return true;
}

bool LayerMRI::SaveVolume( wxWindow* wnd, wxCommandEvent& event )
{
	if ( m_sFilename.size() == 0 || m_volumeRAS.GetPointer() == NULL )
		return false;
		
	m_volumeSource->UpdateMRIFromImage( m_volumeRAS, wnd, event );
	
	wxPostEvent( wnd, event );
	bool bSaved = m_volumeSource->MRIWrite( m_sFilename.c_str() );
	if ( !bSaved )
		m_bModified = true;
	
	return bSaved;
}

void LayerMRI::InitializeVolume()
{
	if ( m_volumeSource == NULL )
		return;

	FSVolume* source = m_volumeSource;
		
	float RASBounds[6];
	source->GetBounds( RASBounds );	
	m_dWorldOrigin[0] = RASBounds[0];
	m_dWorldOrigin[1] = RASBounds[2];
	m_dWorldOrigin[2] = RASBounds[4]; 
	source->GetPixelSize( m_dWorldVoxelSize );
	
	m_dWorldSize[0] = ( ( int )( (RASBounds[1] - RASBounds[0]) / m_dWorldVoxelSize[0] ) ) * m_dWorldVoxelSize[0];
	m_dWorldSize[1] = ( ( int )( (RASBounds[3] - RASBounds[2]) / m_dWorldVoxelSize[1] ) ) * m_dWorldVoxelSize[1];
	m_dWorldSize[2] = ( ( int )( (RASBounds[5] - RASBounds[4]) / m_dWorldVoxelSize[2] ) ) * m_dWorldVoxelSize[2];
	
	m_volumeRAS = source->GetImageOutput();	
}
	
/*
void LayerMRI::InitializeVolume()
{
	if ( m_volumeSource == NULL )
		return;

	FSVolume* source = m_volumeSource;
		
	float RASBounds[6];
	source->GetRASBounds( RASBounds );	
	m_dWorldOrigin[0] = RASBounds[0];
	m_dWorldOrigin[1] = RASBounds[2];
	m_dWorldOrigin[2] = RASBounds[4]; 	
//	MyUtils::GetFSVolumeSourceRASPixelSize( source, m_dWorldVoxelSize );
	source->GetRASPixelSize( m_dWorldVoxelSize );
	
	double sizeX = m_dWorldVoxelSize[0];
	double sizeY = m_dWorldVoxelSize[1];
	double sizeZ = m_dWorldVoxelSize[2];
	
	// Get some values from the MRI.

	m_dWorldSize[0] = RASBounds[1] - RASBounds[0];
	m_dWorldSize[1] = RASBounds[3] - RASBounds[2];
	m_dWorldSize[2] = RASBounds[5] - RASBounds[4];
	
	//
  		// This transforms the voxel space source into RAS space.
	//
	vtkSmartPointer<vtkImageReslice> volumeToRAS = 
			vtkSmartPointer<vtkImageReslice>::New();
	volumeToRAS->SetInput( source->GetOutput() );
	volumeToRAS->SetOutputDimensionality( 3 );

		// This rotates the volume to the proper orientation. From
		// ImageReslice: "applying a transform to the resampling grid (which
		// lies in the output coordinate system) is equivalent to applying
		// the inverse of that transform to the input volume."
//	double rtv[16];
//	source->GetUnscaledRASToVoxelMatrix( rtv );
	double* rtv = source->GetRASToVoxelMatrix();

		// 0,0 0,1 0,2 0,3      rtv[0]  rtv[1]  rtv[2]  0
		// 1,0 1,1 1,2 1,3  =>  rtv[4]  rtv[5]  rtv[6]  0
		// 2,0 2,1 2,2 2,3      rtv[8]  rtv[9]  rtv[10] 0
		// 3,0 3,1 3,2 3,3        0       0       0     1
	
	vtkSmartPointer<vtkMatrix4x4> matrix = 
			vtkSmartPointer<vtkMatrix4x4>::New();
	matrix->SetElement( 0, 0, rtv[0] );
	matrix->SetElement( 0, 1, rtv[1] );
	matrix->SetElement( 0, 2, rtv[2] );
	matrix->SetElement( 0, 3, 0 );
	matrix->SetElement( 1, 0, rtv[4] );
	matrix->SetElement( 1, 1, rtv[5] );
	matrix->SetElement( 1, 2, rtv[6] );
	matrix->SetElement( 1, 3, 0 );
	matrix->SetElement( 2, 0, rtv[8] );
	matrix->SetElement( 2, 1, rtv[9] );
	matrix->SetElement( 2, 2, rtv[10] );
	matrix->SetElement( 2, 3, 0 );
	matrix->SetElement( 3, 0, 0 );
	matrix->SetElement( 3, 1, 0 );
	matrix->SetElement( 3, 2, 0 );
	matrix->SetElement( 3, 3, 1 );

	vtkSmartPointer<vtkTransform> transform = 
			vtkSmartPointer<vtkTransform>::New();
	transform->SetMatrix( matrix );

	volumeToRAS->SetResliceTransform( transform );
	volumeToRAS->BorderOff();

 //  This sets our output extent.
	
	volumeToRAS->SetOutputExtent( 0, ( int )( m_dWorldSize[0] / sizeX - 1 ), 
								  0, ( int )( m_dWorldSize[1] / sizeY - 1 ),
								  0, ( int )( m_dWorldSize[2] / sizeZ - 1 ) );
	
	volumeToRAS->SetOutputSpacing( sizeX, sizeY, sizeZ );
//	volumeToRAS->SetOutputOrigin( RASBounds[0], RASBounds[2], RASBounds[4] );
	
	m_dWorldSize[0] = ( ( int )( m_dWorldSize[0] / sizeX ) ) * sizeX;
	m_dWorldSize[1] = ( ( int )( m_dWorldSize[1] / sizeY ) ) * sizeY;
	m_dWorldSize[2] = ( ( int )( m_dWorldSize[2] / sizeZ ) ) * sizeZ;
	
//	cout << "new world size: " << mWorldSize[0] << " " << mWorldSize[1] << " " << mWorldSize[2] << endl;
		
	volumeToRAS->Update();
	source->SetOriginalOrigin( volumeToRAS->GetOutput()->GetOrigin() );
	
	vtkSmartPointer<vtkImageChangeInformation> infoFilter = vtkSmartPointer<vtkImageChangeInformation>::New();
	infoFilter->SetInput( volumeToRAS->GetOutput() );
	infoFilter->SetOutputOrigin( RASBounds[0], RASBounds[2], RASBounds[4] );
	m_volumeRAS = infoFilter->GetOutput();	
}
*/		
		
void LayerMRI::InitializeActors()
{
	for ( int i = 0; i < 3; i++ )
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
	this->UpdateColorMap();
}

void LayerMRI::UpdateOpacity()
{
	for ( int i = 0; i < 3; i++ )
	{
		m_sliceActor2D[i]->SetOpacity( mProperties->GetOpacity() );
		m_sliceActor3D[i]->SetOpacity( mProperties->GetOpacity() );
	}	
}

void LayerMRI::UpdateColorMap () 
{
	assert( mProperties );

	for ( int i = 0; i < 3; i++ )
		mColorMap[i]->SetActiveComponent( m_nActiveFrame );
	
	switch ( mProperties->GetColorMap() ) 
	{
		case LayerPropertiesMRI::NoColorMap:
			for ( int i = 0; i < 3; i++ )
				mColorMap[i]->SetLookupTable( NULL );
			break;
    
		case LayerPropertiesMRI::Grayscale:
			for ( int i = 0; i < 3; i++ )
				mColorMap[i]->SetLookupTable( mProperties->GetGrayScaleTable() );
			break;
    
		case LayerPropertiesMRI::Heat:
			for ( int i = 0; i < 3; i++ )
					mColorMap[i]->SetLookupTable( mProperties->GetHeatScaleTable() );
			break;
		case LayerPropertiesMRI::Jet:
			for ( int i = 0; i < 3; i++ )
				mColorMap[i]->SetLookupTable( mProperties->GetJetScaleTable() );
			break;
		case LayerPropertiesMRI::LUT:
			for ( int i = 0; i < 3; i++ )
				mColorMap[i]->SetLookupTable( mProperties->GetLUTTable() );
			break;
		default:
			break;
	}
}

void LayerMRI::UpdateResliceInterpolation () 
{
	assert( mProperties );

	for ( int i = 0; i < 3; i++ )
	{
		if( mReslice[i].GetPointer() ) 
		{
			mReslice[i]->SetInterpolationMode( mProperties->GetResliceInterpolation() );
		}
	}
}

void LayerMRI::UpdateTextureSmoothing () 
{
	assert( mProperties );

	for ( int i = 0; i < 3; i++ )
	{
		m_sliceActor2D[i]->SetInterpolate( mProperties->GetTextureSmoothing() );
		m_sliceActor3D[i]->SetInterpolate( mProperties->GetTextureSmoothing() );
	}
}

void LayerMRI::Append2DProps( vtkRenderer* renderer, int nPlane )
{
	wxASSERT ( nPlane >= 0 && nPlane <= 2 );
	
	renderer->AddViewProp( m_sliceActor2D[nPlane] );
}

void LayerMRI::Append3DProps( vtkRenderer* renderer )
{
	for ( int i = 0; i < 3; i++ )
		renderer->AddViewProp( m_sliceActor3D[i] ); 
}

/*
void LayerMRI::SetSliceNumber( int* sliceNumber )
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

void LayerMRI::SetSlicePositionToWorldCenter()
{
	if ( m_volumeSource == NULL )
		return;
	
	// Get some values from the MRI.
	double pos[3];
	for ( int i = 0; i < 3; i++ )
		pos[i] = ((int)( m_dWorldSize[i]/2/m_dWorldVoxelSize[i] ) + 0.0 ) * m_dWorldVoxelSize[i] + m_dWorldOrigin[i];
	
	SetSlicePosition( pos );
}

void LayerMRI::OnSlicePositionChanged( int nPlane ) 
{  
	if ( m_volumeSource == NULL )
		return;
	
	assert( mProperties );
	
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
}

LayerPropertiesMRI* LayerMRI::GetProperties()
{
	return mProperties;
}

void LayerMRI::DoListenToMessage( std::string const iMessage, void* const iData )
{
	if ( iMessage == "ColorMapChanged" )
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
	}
}

void LayerMRI::SetVisible( bool bVisible )
{
	for ( int i = 0; i < 3; i++ )
	{
		m_sliceActor2D[i]->SetVisibility( bVisible ? 1 : 0 );
		m_sliceActor3D[i]->SetVisibility( bVisible ? 1 : 0 );
	}
	this->SendBroadcast( "LayerActorUpdated", this );
}

bool LayerMRI::IsVisible()
{
	return m_sliceActor2D[0]->GetVisibility() > 0;
}

double LayerMRI::GetVoxelValue( double* pos )
{
	if ( m_volumeRAS.GetPointer() == NULL )
		return 0;
	
	double* orig = m_volumeRAS->GetOrigin();
	double* vsize = m_volumeRAS->GetSpacing();
	int* ext = m_volumeRAS->GetExtent();
	
	int n[3];
	for ( int i = 0; i < 3; i++ )
	{
		n[i] = (int)( ( pos[i] - orig[i] ) / vsize[i] + 0.5 );
	}
	
	if ( n[0] < ext[0] || n[0] > ext[1] || 
		 n[1] < ext[2] || n[1] > ext[3] || 
		 n[2] < ext[4] || n[2] > ext[5] )
		return 0;
	else
		return m_volumeRAS->GetScalarComponentAsDouble( n[0], n[1], n[2], m_nActiveFrame );
}

void LayerMRI::SetModified()
{
	mReslice[0]->Modified();
	mReslice[1]->Modified();
	mReslice[2]->Modified();
	
	LayerEditable::SetModified();
}

std::string LayerMRI::GetLabelName( double value )
{
	int nIndex = (int)value;
	if ( GetProperties()->GetColorMap() == LayerPropertiesMRI::LUT )
	{
		COLOR_TABLE* ct = GetProperties()->GetLUTCTAB();
		if ( !ct )
			return "";
		char name[128];
		int nValid = 0;
		int nTotalCount = 0;
		CTABgetNumberOfTotalEntries( ct, &nTotalCount );
		if ( nIndex < nTotalCount )
			CTABisEntryValid( ct, nIndex, &nValid );
		if ( nValid && CTABcopyName( ct, nIndex, name, 128 ) == 0 )
		{
			return name;
		}	
	}
	
	return "";
}

void LayerMRI::RemapPositionToRealRAS( const double* pos_in, double* pos_out )
{
	m_volumeSource->RemapPositionToRealRAS( pos_in, pos_out );
}

void LayerMRI::RemapPositionToRealRAS( double x_in, double y_in, double z_in, 
							double& x_out, double& y_out, double& z_out )
{
	m_volumeSource->RemapPositionToRealRAS( x_in, y_in, z_in, x_out, y_out, z_out );
}

int LayerMRI::GetNumberOfFrames()
{
	if ( m_volumeRAS )
		return m_volumeRAS->GetNumberOfScalarComponents();
	else
		return 1;
}

void LayerMRI::SetActiveFrame( int nFrame )
{
	if ( nFrame != m_nActiveFrame && nFrame >= 0 && nFrame < this->GetNumberOfFrames() )
	{
		m_nActiveFrame = nFrame;
		this->DoListenToMessage( "ColorMapChanged", this );		
		this->SendBroadcast( "LayerActiveFrameChanged", this );
	}
}

