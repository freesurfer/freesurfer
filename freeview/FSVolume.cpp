/**
 * @file  FSVolume.h
 * @brief Interactor to manage mouse and key input on render view.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2008/06/04 20:43:24 $
 *    $Revision: 1.3.2.1 $
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
#include "FSVolume.h"
#include <stdexcept>
#include "vtkShortArray.h"
#include "vtkLongArray.h"
#include "vtkUnsignedCharArray.h"
#include "vtkFloatArray.h"
#include "vtkPointData.h"
#include "vtkCellData.h"
#include "vtkIntArray.h"
#include "vtkSmartPointer.h"
#include "vtkImageReslice.h"
#include "vtkMatrix4x4.h"
#include "vtkTransform.h"
#include "vtkImageChangeInformation.h"

using namespace std;

FSVolume::FSVolume() :
	m_MRI( NULL ),
	m_fMinValue( 0 ),
	m_fMaxValue( 1 ),
	m_bResampleToRAS( true ),
	m_bBoundsCacheDirty( true )
{
}
	
FSVolume::~FSVolume()
{
	if ( m_MRI )
		::MRIfree( &m_MRI );
}
		
bool FSVolume::MRIRead( const char* filename, wxWindow* wnd, wxCommandEvent& event  )
{
	if ( m_MRI )
		::MRIfree( &m_MRI );
	
	event.SetInt( event.GetInt() + 1 );
	wxPostEvent( wnd, event );
	
	char* fn = strdup( filename );
	m_MRI = ::MRIread( fn );
	free( fn );	
	
	if ( m_MRI == NULL ) {
		cerr << "MRIread failed";
		return false;
	}	
	
	MRIvalRange( m_MRI, &m_fMinValue, &m_fMaxValue );
	
	this->CopyMatricesFromMRI();
	this->MapMRIToImage( wnd, event );
	
	return true;
}

void FSVolume::Create( FSVolume* src_vol, bool bCopyVoxelData )
{
	if ( m_MRI )
		::MRIfree( &m_MRI );
	
	memcpy( m_RASToVoxelMatrix, src_vol->m_RASToVoxelMatrix, 16 * sizeof( double ) );
	memcpy( m_VoxelToRASMatrix, src_vol->m_VoxelToRASMatrix, 16 * sizeof( double ) );
	memcpy( m_MRIToImageMatrix, src_vol->m_MRIToImageMatrix, 16 * sizeof( double ) );
	
	m_fMinValue = src_vol->m_fMinValue;
	m_fMaxValue = src_vol->m_fMaxValue;
	m_bResampleToRAS = src_vol->m_bResampleToRAS; 
	
	SetOriginalOrigin( src_vol->m_fOriginalOrigin );
	
	if ( bCopyVoxelData )
	{
		m_MRI = MRIcopy( src_vol->m_MRI, NULL );
	}
	else
	{
		m_MRI = MRIallocSequence( src_vol->m_MRI->width, src_vol->m_MRI->height, src_vol->m_MRI->depth,
							  MRI_FLOAT, src_vol->m_MRI->nframes );
	}
	
	if ( NULL == m_MRI ) 
	{
		cerr << "Couldn't allocate new mri." << endl;
		return;
	}

  	// Copy the header from the template into the new mri.
	if ( !bCopyVoxelData )
		MRIcopyHeader( src_vol->m_MRI, m_MRI );
	
	m_imageData = src_vol->m_imageData;
}
	
bool FSVolume::MRIWrite( const char* filename )
{
	char* fn = strdup( filename );
	int err = ::MRIwrite( m_MRI, fn );
	free( fn );

	if ( err != 0 ) 
	{
		cerr << "MRIwrite failed" << endl;
	}
	
	return err == 0;
}

bool FSVolume::MRIWrite()
{
	int err = ::MRIwrite( m_MRI, m_MRI->fname );
	if ( err != 0 ) 
	{
		cerr << "MRIwrite failed" << endl;
	}
	return err == 0;
}

void FSVolume::UpdateMRIFromImage( vtkImageData* rasImage, wxWindow* wnd, wxCommandEvent& event )				
{
	vtkSmartPointer<vtkImageReslice> rasToNative = 
			vtkSmartPointer<vtkImageReslice>::New();
	
	vtkSmartPointer<vtkImageChangeInformation> infoFilter = 
			vtkSmartPointer<vtkImageChangeInformation>::New();
	
	infoFilter->SetInput( rasImage );
	infoFilter->SetOutputOrigin( m_fOriginalOrigin );
	rasToNative->SetInput( infoFilter->GetOutput() );
	rasToNative->SetOutputDimensionality( 3 );

/*	double* rtv = GetVoxelToRASMatrix();
	vtkSmartPointer<vtkMatrix4x4> matrix = vtkSmartPointer<vtkMatrix4x4>::New();
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
	matrix->SetElement( 3, 3, 1 );*/
	double inv[16];
	vtkMatrix4x4::Invert( m_MRIToImageMatrix, inv );
	vtkSmartPointer<vtkMatrix4x4> matrix = vtkSmartPointer<vtkMatrix4x4>::New();
	matrix->DeepCopy( inv );

	vtkSmartPointer<vtkTransform> transform = 
			vtkSmartPointer<vtkTransform>::New();
	transform->SetMatrix( matrix );

	rasToNative->SetResliceTransform( transform );
	rasToNative->BorderOff();
	rasToNative->SetOutputExtent( 0, m_MRI->width - 1, 
								  0, m_MRI->height - 1,
								  0, m_MRI->depth - 1 );

	if ( m_bResampleToRAS )
		rasToNative->SetOutputSpacing( 1, 1, 1 );
	else
		rasToNative->SetOutputSpacing( m_MRI->xsize, m_MRI->ysize, m_MRI->zsize );
	rasToNative->Update();
	vtkImageData* out = rasToNative->GetOutput();
	
	int nProgressStep = ( 99 - event.GetInt() ) / 5;
	for ( int j = 0; j < m_MRI->height; j++ )
	{
		for ( int k = 0; k < m_MRI->depth; k++ )
		{
			void* ptr = out->GetScalarPointer( 0, j, k );
			BUFTYPE* buf = &MRIseq_vox( m_MRI, 0, j, k, 0);
			memcpy( buf, ptr, m_MRI->width*out->GetScalarSize() );		
		}
		if ( m_MRI->height >= 5 && j%(m_MRI->height/5) == 0 )
		{
			event.SetInt( event.GetInt() + nProgressStep );
			wxPostEvent( wnd, event );
		}
	}
	MRIvalRange( m_MRI, &m_fMinValue, &m_fMaxValue );
}
		
int FSVolume::OriginalIndexToRAS( float iIdxX, float iIdxY, float iIdxZ,
					float& oRASX, float& oRASY, float& oRASZ )
{
	if ( m_MRI == NULL ) {
		cerr << "No MRI is present." << endl;
		return 1;
	}

	Real ix, iy, iz, wx, wy, wz;
	int r;

	ix = iIdxX;
	iy = iIdxY;
	iz = iIdxZ;
	r = ::MRIvoxelToWorld( m_MRI, ix, iy, iz, &wx, &wy, &wz );
	oRASX = wx;
	oRASY = wy;
	oRASZ = wz;

	return r;
}

int FSVolume::RASToOriginalIndex ( float iRASX, float iRASY, float iRASZ,
					float& oIdxX, float& oIdxY, float& oIdxZ )
{
	if ( m_MRI == NULL ) {
		cerr << "No MRI is present." << endl;
		return 1;
	}
	
	Real ix, iy, iz, wx, wy, wz;
	int r;

	wx = iRASX;
	wy = iRASY;
	wz = iRASZ;
	r = ::MRIworldToVoxel( m_MRI, wx, wy, wz, &ix, &iy, &iz );
	oIdxX = ix;
	oIdxY = iy;
	oIdxZ = iz;

	return r;
}
	
int FSVolume::RASToOriginalIndex ( float iRASX, float iRASY, float iRASZ,
					int& oIdxX, int& oIdxY, int& oIdxZ )
{
	if ( m_MRI == NULL ) {
		cerr << "No MRI is present." << endl;
		return 1;
	}
	
	Real wx, wy, wz;
	int ix, iy, iz;
	int r;

	wx = iRASX;
	wy = iRASY;
	wz = iRASZ;
	r = ::MRIworldToVoxelIndex( m_MRI, wx, wy, wz, &ix, &iy, &iz );
	oIdxX = ix;
	oIdxY = iy;
	oIdxZ = iz;

	return r;
}


void FSVolume::MapMRIToImage ( wxWindow* wnd, wxCommandEvent& event ) 
{
	event.SetInt( event.GetInt() + 1 );
	wxPostEvent( wnd, event );
	
	// first copy mri data to image
	vtkDataArray *scalars = NULL;
	vtkUnsignedCharArray  *ucharScalars = NULL;
	vtkIntArray           *intScalars = NULL;
	vtkShortArray         *shortScalars = NULL;
	vtkLongArray          *longScalars = NULL;
	vtkFloatArray         *floatScalars = NULL;
	void* pixels;
	int cValues;
	int nZ, nY;
	int cRead;
	int zElement=0;

	if ( m_MRI == NULL ) 
	{
		cerr << "No MRI is present." << endl;
		return;
	}

	int zX = m_MRI->width;
	int zY = m_MRI->height;
	int zZ = m_MRI->depth;
	int zFrames = m_MRI->nframes;

	vtkImageData* imageData = vtkImageData::New();
	
  // This object's output space is in voxel coordinates.
	imageData->SetDimensions( zX, zY, zZ );
  
	// if map to RAS
	if ( m_bResampleToRAS )
	{
		imageData->SetSpacing( 1, 1, 1 );
		imageData->SetOrigin( -zX/2, -zY/2, -zZ/2 );
	//	imageData->SetOrigin( 0, 0, 0 );
	}
	else
	{ 
		imageData->SetSpacing( m_MRI->xsize, m_MRI->ysize, m_MRI->zsize );
		imageData->SetOrigin( -zX/2*m_MRI->xsize, -zY/2*m_MRI->ysize, -zZ/2*m_MRI->zsize );
	}  
	
	imageData->SetWholeExtent( 0, zX-1, 0, zY-1, 0, zZ-1 );
	imageData->SetNumberOfScalarComponents( zFrames );

  // create the scalars for all of the images. set the element size
  // for the data we will read.
	switch ( m_MRI->type ) {
		case MRI_UCHAR:
			imageData->SetScalarTypeToUnsignedChar();
			ucharScalars = vtkUnsignedCharArray::New();
			scalars = (vtkDataArray*) ucharScalars;
			zElement = sizeof( unsigned char );
			break;
		case MRI_INT:
			imageData->SetScalarTypeToInt();
			intScalars = vtkIntArray::New();
			scalars = (vtkDataArray*) intScalars;
			zElement = sizeof( int );
			break;
		case MRI_LONG:
			imageData->SetScalarTypeToLong();
			longScalars = vtkLongArray::New();
			scalars = (vtkDataArray*) longScalars;
			zElement = sizeof( long );
			break;
		case MRI_FLOAT:
			imageData->SetScalarTypeToFloat();
			floatScalars = vtkFloatArray::New();
			scalars = (vtkDataArray*) floatScalars;
			zElement = sizeof( float );
			break;
		case MRI_SHORT:
			imageData->SetScalarTypeToShort();
			shortScalars = vtkShortArray::New();
			scalars = (vtkDataArray*) shortScalars;
			zElement = sizeof( short );
			break;
		default:
			break ;
	}

	if ( NULL == scalars ) 
	{
		cerr << "Couldn't allocate scalars array." << endl;
		return;
	}
  
  // change the number of components to store tuples
	if( zFrames > 1 ) 
	{
		scalars->SetNumberOfComponents( zFrames );
	}

	cValues = zX * zY * zZ;
	scalars->Allocate( cValues );

	event.SetInt( event.GetInt() + ( 100 - event.GetInt() ) / 10 );
	wxPostEvent( wnd, event );
	int nProgressStep = ( 99-event.GetInt() ) / 5;
	
  // Copy the slice data into the scalars.
//	vtkDebugMacro (<< "Copying " << cValues << " values into the scalars array");
//	float* tuple = (float*) malloc( sizeof(float) * zFrames );
	cRead = 0;
	int nTuple = 0;
	pixels = NULL;
	scalars->SetNumberOfTuples( zX*zY*zZ );
	for ( nZ = 0; nZ < zZ; nZ++ ) 
	{
		for ( nY = 0; nY < zY; nY++ ) {
#if 1
      for ( int nX = 0; nX < zX; nX++ ) {
	for ( int nFrame = 0; nFrame < zFrames; nFrame++ ) {
	  
		switch ( m_MRI->type ) {
			case MRI_UCHAR:
				scalars->SetComponent( nTuple, nFrame, 
										  MRIseq_vox( m_MRI, nX, nY, nZ, nFrame ) );
				break;
			case MRI_INT:
				scalars->SetComponent( nTuple, nFrame, 
										  MRIIseq_vox( m_MRI, nX, nY, nZ, nFrame ) );
				break;
			case MRI_LONG:
				scalars->SetComponent( nTuple, nFrame, 
										  MRILseq_vox( m_MRI, nX, nY, nZ, nFrame ) );
				break;
			case MRI_FLOAT:
				scalars->SetComponent( nTuple, nFrame, 
										  MRIFseq_vox( m_MRI, nX, nY, nZ, nFrame ) );
				break;
			case MRI_SHORT:
				scalars->SetComponent( nTuple, nFrame, 
										  MRISseq_vox( m_MRI, nX, nY, nZ, nFrame ) );                    
				break;
			default:
				break;
		}
	  
	}
	nTuple++;
                
	  }
      
#else
      // Not sure how this will work when reading in frames.
//      vtkDebugMacro (<< "Getting a write pointer for " <<
//                     cValues << " values");
      switch ( m_MRI->type ) {
		  case MRI_UCHAR:
			  pixels = (void*) ucharScalars->WritePointer( cRead, zX * zElement );
			  break;
		  case MRI_INT:
			  pixels = (void*) intScalars->WritePointer( cRead, zX * zElement );
			  break;
		  case MRI_LONG:
			  pixels = (void*) longScalars->WritePointer( cRead, zX * zElement );
			  break;
		  case MRI_FLOAT:
			  pixels = (void*) floatScalars->WritePointer( cRead, zX * zElement );
			  break;
		  case MRI_SHORT:
			  pixels = (void*) shortScalars->WritePointer( cRead, zX * zElement );
			  break;
		  default:
			  break ;
	  }
	  if ( NULL == pixels ) {
		  cerr << "Couldn't get a write pointer" << endl;
		  return;
	  }
	  memcpy( pixels, m_MRI->slices[nZ][nY], zX* zElement );
	  cRead += zX;
#endif
		}
	
		if ( nZ%(max(1, zZ/5)) == 0 )
		{
			event.SetInt( event.GetInt() + nProgressStep );  
			wxPostEvent( wnd, event );
		}
	}
//	free( tuple );
  
  // Assign the scalars array to the image.
	imageData->GetPointData()->SetScalars( scalars );  
	scalars->Delete();

  // Need to recalc our bounds at some point.
	m_bBoundsCacheDirty = true;
	

	/////////////////////////////////////////////////////
	// then remap to ras space

	vtkSmartPointer<vtkImageReslice> volumeToRAS = 
			vtkSmartPointer<vtkImageReslice>::New();
	volumeToRAS->SetInput( imageData );
	volumeToRAS->SetOutputDimensionality( 3 );
	imageData->Delete();
	
	double* rtv = this->GetRASToVoxelMatrix();
//	double* rtv = this->GetVoxelToRASMatrix();
	
		// 0,0 0,1 0,2 0,3      rtv[0]  rtv[1]  rtv[2]  0
		// 1,0 1,1 1,2 1,3  =>  rtv[4]  rtv[5]  rtv[6]  0
		// 2,0 2,1 2,2 2,3      rtv[8]  rtv[9]  rtv[10] 0
		// 3,0 3,1 3,2 3,3        0       0       0     1

	vtkSmartPointer<vtkMatrix4x4> matrix = 
			vtkSmartPointer<vtkMatrix4x4>::New();
	matrix->Zero();
	vtkSmartPointer<vtkTransform> transform = 
			vtkSmartPointer<vtkTransform>::New();

	volumeToRAS->SetResliceTransform( transform );
	volumeToRAS->BorderOff();
	
	vtkSmartPointer<vtkImageChangeInformation> infoFilter = vtkSmartPointer<vtkImageChangeInformation>::New();	
	double rasVoxelSize[3];
	this->GetPixelSize( rasVoxelSize );
	volumeToRAS->SetOutputSpacing( rasVoxelSize );	
	if ( m_bResampleToRAS )
	{
		float RASBounds[6];
		this->GetBounds( RASBounds );
		
		matrix->SetElement( 0, 0, rtv[0] );
		matrix->SetElement( 0, 1, rtv[1] );
		matrix->SetElement( 0, 2, rtv[2] );
		
		matrix->SetElement( 1, 0, rtv[4] );
		matrix->SetElement( 1, 1, rtv[5] );
		matrix->SetElement( 1, 2, rtv[6] );
		
		matrix->SetElement( 2, 0, rtv[8] );
		matrix->SetElement( 2, 1, rtv[9] );
		matrix->SetElement( 2, 2, rtv[10] );
		
		matrix->SetElement( 3, 3, 1 );

	//	matrix->DeepCopy( rtv );
		
		transform->SetMatrix( matrix );
		volumeToRAS->SetOutputExtent( 0, ( int )( (RASBounds[1] - RASBounds[0]) / rasVoxelSize[0] - 1 ), 
									0, ( int )( (RASBounds[3] - RASBounds[2]) / rasVoxelSize[1] - 1 ),
									0, ( int )( (RASBounds[5] - RASBounds[4]) / rasVoxelSize[2] - 1 ) );	
	
		infoFilter->SetOutputOrigin( RASBounds[0], RASBounds[2], RASBounds[4] );
		
//		cout << "Origin_2: " << RASBounds[0] << " " << RASBounds[2] << " " << RASBounds[4] << endl;
	}
	else
	{
		int nExt[3] = { 1 };
		if ( fabs( rtv[0] ) > fabs( rtv[1] ) && fabs( rtv[0] ) > fabs( rtv[2] ) )
		{
			matrix->SetElement( 0, 0, rtv[0] > 0 ? 1 : -1 );
			nExt[0] = zX;
		}
		else if ( fabs( rtv[1] ) > fabs( rtv[2] ) )
		{
			matrix->SetElement( 0, 1, rtv[1] > 0 ? 1 : -1 );
			nExt[1] = zX;
		}
		else 
		{
			matrix->SetElement( 0, 2, rtv[2] > 0 ? 1 : -1 );
			nExt[2] = zX;
		}
		
		if ( fabs( rtv[4] ) > fabs( rtv[5] ) && fabs( rtv[4] ) > fabs( rtv[6] ) )
		{
			matrix->SetElement( 1, 0, rtv[4] > 0 ? 1 : -1 );
			nExt[0] = zY;
		}
		else if ( fabs( rtv[5] ) > fabs( rtv[6] ) )
		{
			matrix->SetElement( 1, 1, rtv[5] > 0 ? 1 : -1 );
			nExt[1] = zY;
		}
		else 
		{
			matrix->SetElement( 1, 2, rtv[6] > 0 ? 1 : -1 );
			nExt[2] = zY;
		}
		
		if ( fabs( rtv[8] ) > fabs( rtv[9] ) && fabs( rtv[8] ) > fabs( rtv[10] ) )
		{
			matrix->SetElement( 2, 0, rtv[8] > 0 ? 1 : -1 );
			nExt[0] = zZ;
		}
		else if ( fabs( rtv[9] ) > fabs( rtv[10] ) )
		{
			matrix->SetElement( 2, 1, rtv[9] > 0 ? 1 : -1 );
			nExt[1] = zZ;
		}
		else 
		{
			matrix->SetElement( 2, 2, rtv[10] > 0 ? 1 : -1 );
			nExt[2] = zZ;
		}
		
		matrix->SetElement( 3, 3, 1 );
		
		transform->SetMatrix( matrix );
		volumeToRAS->SetOutputExtent( 0, nExt[0]-1, 
									  0, nExt[1]-1, 
									  0, nExt[2]-1 );	
		infoFilter->SetOutputOrigin( 0, 0, 0 );
	}
	infoFilter->SetInput( volumeToRAS->GetOutput() );
	infoFilter->Update();
	SetOriginalOrigin( volumeToRAS->GetOutput()->GetOrigin() );		
	m_imageData = infoFilter->GetOutput();	
	vtkMatrix4x4::DeepCopy( m_MRIToImageMatrix, matrix );
}


vtkImageData* FSVolume::GetImageOutput()
{
	return m_imageData;
}

void FSVolume::CopyMatricesFromMRI () 
{
	if ( NULL == m_MRI )
		return;

	MATRIX* m = extract_i_to_r( m_MRI );
	for ( int i = 0; i < 16; i++ ) 
	{
		m_VoxelToRASMatrix[i] = (double) *MATRIX_RELT((m),(i/4)+1,(i%4)+1);
	}
	MatrixFree( &m );

	m = extract_r_to_i( m_MRI );
	for ( int i = 0; i < 16; i++ ) 
	{
		m_RASToVoxelMatrix[i] = (double) *MATRIX_RELT((m),(i/4)+1,(i%4)+1);
	}
	MatrixFree( &m );
}

void FSVolume::GetBounds ( float oRASBounds[6] ) 
{
	if ( NULL == m_MRI ) {

		oRASBounds[0] = oRASBounds[1] = oRASBounds[2] =
				oRASBounds[3] = oRASBounds[4] = oRASBounds[5] = 0;

		return;
	}
	
	if ( !m_bResampleToRAS )
	{
		double* orig = m_imageData->GetOrigin();
		double* spacing = m_imageData->GetSpacing();
		int* dim = m_imageData->GetDimensions();
		for ( int i = 0; i < 3; i++ )
		{
			oRASBounds[i*2] = ( float )orig[i];
			oRASBounds[i*2+1] = ( float )( orig[i] + spacing[i] * dim[i] );
		} 
		return;
	}

  // If we need to rebuild the cache...
	if ( m_bBoundsCacheDirty ) 
	{
		m_RASBounds[0] = m_RASBounds[2] = m_RASBounds[4] = 999999;
		m_RASBounds[1] = m_RASBounds[3] = m_RASBounds[5] = -999999;

    // For each corner, convert to RAS, and get the bounds.
		float cornerFactor[3];
		float RAS[3];
		for ( cornerFactor[2] = 0; cornerFactor[2] <= 1; cornerFactor[2]++ ) {
			for ( cornerFactor[1] = 0; cornerFactor[1] <= 1; cornerFactor[1]++ ) {
				for ( cornerFactor[0] = 0; cornerFactor[0] <= 1; cornerFactor[0]++ ) {

					this->OriginalIndexToRAS( cornerFactor[0]*m_MRI->width,
											 cornerFactor[1]*m_MRI->height,
											 cornerFactor[2]*m_MRI->depth,
											 RAS[0], RAS[1], RAS[2] );

					if ( RAS[0] < m_RASBounds[0] ) m_RASBounds[0] = RAS[0];
					if ( RAS[0] > m_RASBounds[1] ) m_RASBounds[1] = RAS[0];
					if ( RAS[1] < m_RASBounds[2] ) m_RASBounds[2] = RAS[1];
					if ( RAS[1] > m_RASBounds[3] ) m_RASBounds[3] = RAS[1];
					if ( RAS[2] < m_RASBounds[4] ) m_RASBounds[4] = RAS[2];
					if ( RAS[2] > m_RASBounds[5] ) m_RASBounds[5] = RAS[2];
				}
			}
		}

		m_bBoundsCacheDirty = false;
	}

	oRASBounds[0] = m_RASBounds[0];
	oRASBounds[1] = m_RASBounds[1];
	oRASBounds[2] = m_RASBounds[2];
	oRASBounds[3] = m_RASBounds[3];
	oRASBounds[4] = m_RASBounds[4];
	oRASBounds[5] = m_RASBounds[5];
}

int FSVolume::GetNumberOfFrames()
{
	if ( m_MRI == NULL ) 
	{
		cerr << "No MRI is present." << endl;
		return 0;
	}
	
	return m_MRI->nframes;
}

void FSVolume::GetPixelSize( double* pixelSize )
{
	if ( m_MRI == NULL ) 
	{
		cerr << "No MRI is present." << endl;
		return;
	}
		
	double* m = GetVoxelToRASMatrix();
		
	if ( fabs( m[0] ) > fabs( m[1] ) && fabs( m[0] ) > fabs( m[2] ) )
		pixelSize[0] = m_MRI->xsize;
	else if ( fabs( m[1] ) > fabs( m[2] ) )
		pixelSize[0] = m_MRI->ysize;
	else
		pixelSize[0] = m_MRI->zsize;
		
	if ( fabs( m[4] ) > fabs( m[5] ) && fabs( m[4] ) > fabs( m[6] ) )
		pixelSize[1] = m_MRI->xsize;
	else if ( fabs( m[5] ) > fabs( m[6] ) )
		pixelSize[1] = m_MRI->ysize;
	else
		pixelSize[1] = m_MRI->zsize;
		
	if ( fabs( m[8] ) > fabs( m[9] ) && fabs( m[8] ) > fabs( m[10] ) )
		pixelSize[2] = m_MRI->xsize;
	else if ( fabs( m[9] ) > fabs( m[10] ) )
		pixelSize[2] = m_MRI->ysize;
	else
		pixelSize[2] = m_MRI->zsize;
}

void FSVolume::SetOriginalOrigin( double* origin )
{
	for ( int i = 0; i < 3; i++ )
		m_fOriginalOrigin[i] = origin[i];
}

void FSVolume::RemapPositionToRealRAS( double x_in, double y_in, double z_in, 
							 double& x_out, double& y_out, double& z_out )
{
	double pos[3] = { x_in, y_in, z_in };
	RemapPositionToRealRAS( pos, pos );
	x_out = pos[0];
	y_out = pos[1];
	z_out = pos[2];
}

void FSVolume::RemapPositionToRealRAS( const double* pos_in, double* pos_out )
{
	// if already resampled to standard RAS, no need to remap
	if ( m_bResampleToRAS )
		memcpy( pos_out, pos_in, sizeof( double ) * 3 );
	else
	{
		double pos[4] = { 0 };	
		double ident[4] = { 1, 1, 1, 0 };	
		double* vs = m_imageData->GetSpacing();
		for ( int i = 0; i < 3; i++ )
			pos[i] = pos_in[i] / vs[i];
//		double mat[16];
//		vtkMatrix4x4::Transpose( m_MRIToImageMatrix, mat );
		vtkMatrix4x4::MultiplyPoint(m_MRIToImageMatrix, pos, pos );
		vtkMatrix4x4::MultiplyPoint(m_MRIToImageMatrix, ident, ident );
//		cout << "remap: " << pos[0] << " " << pos[1] << " " << pos[2] << endl;
		int dim[3] = { m_MRI->width, m_MRI->height, m_MRI->depth };
		for ( int i = 0; i < 3; i++ )
		{
			if ( ident[i] < 0 )
				pos[i] = dim[i] + pos[i];
		}
		float fpos[3];
		OriginalIndexToRAS( (float)pos[0], (float)pos[1], (float)pos[2], fpos[0], fpos[1], fpos[2] );
//		cout << "out: " << fpos[0] << " " << fpos[1] << " " << fpos[2] << endl;
		for ( int i = 0; i < 3; i++ )
		{
			pos_out[i] = fpos[i];
		}
	}
}

void FSVolume::RASToOutputIndex( const double* pos_in, int* index_out )
{
	if ( m_bResampleToRAS )
	{
		double* vs = m_imageData->GetSpacing();
		double* orig = m_imageData->GetOrigin();
		for ( int i = 0; i < 3; i++ )
		{
			index_out[i] = ( int )( ( pos_in[i] - orig[i] ) / vs[i] + 0.5 );
		}
	//	cout << vs[0] << " " << vs[1] << " " << vs[2] << endl;
	}
	else
	{
		float pos[4] = { 0 };
		RASToOriginalIndex( pos_in[0], pos_in[1], pos_in[2], pos[0], pos[1], pos[2] );
		double ident[4] = { 1, 1, 1, 0 };
		double mat[16];
		vtkMatrix4x4::Transpose( m_MRIToImageMatrix, mat );
		vtkMatrix4x4::MultiplyPoint(mat, pos, pos );
		vtkMatrix4x4::MultiplyPoint(mat, ident, ident );
//		cout << "remap: " << pos[0] << " " << pos[1] << " " << pos[2] << endl;
		int* dim = m_imageData->GetDimensions();
		for ( int i = 0; i < 3; i++ )
		{
			if ( ident[i] < 0 )
				pos[i] = dim[i] + pos[i];
			
			index_out[i] = ( int )( pos[i] + 0.5 );
		}
	}
}
