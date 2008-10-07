/**
 * @file  LayerDTI.cpp
 * @brief Layer data object for MRI volume.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2008/10/07 22:01:54 $
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
#include "LayerDTI.h"
#include "LayerPropertiesDTI.h"
#include "MyUtils.h"
#include "FSVolume.h"
#include "vtkFloatArray.h"
#include "vtkPointData.h"
#include "vtkImageMapToColors.h"
#include "vtkLookupTable.h"

LayerDTI::LayerDTI( LayerMRI* ref ) : LayerMRI( ref ),
		m_vectorSource( NULL)
{
	m_strTypeNames.push_back( "DTI" );
	if ( mProperties )
		delete mProperties;
	
	mProperties = new LayerPropertiesDTI();
	mProperties->AddListener( this );
	
	SetEditable( false );
}

LayerDTI::~LayerDTI()
{	
	if ( m_vectorSource )
		delete m_vectorSource;
}

LayerPropertiesDTI*	LayerDTI::GetProperties()
{
	return ( LayerPropertiesDTI* )mProperties;
}

bool LayerDTI::LoadDTIFromFile( wxWindow* wnd, wxCommandEvent& event )
{
	if ( !LayerMRI::LoadVolumeFromFile( wnd, event ) )
		return false;
	
	if ( m_vectorSource )
		delete m_vectorSource;
	
	m_vectorSource = new FSVolume( m_volumeRef );
	m_vectorSource->SetResampleToRAS( m_bResampleToRAS );
	event.SetInt( 25 );
	if ( !m_vectorSource->MRIRead( 	m_sVectorFileName.c_str(),  
		  							m_sRegFilename.size() > 0 ? m_sRegFilename.c_str() : NULL,
									wnd, 
									event ) )
		return false;
	
	if ( m_vectorSource->GetNumberOfFrames() < 3 )
	{
		SetErrorString( "Vector data file is not valid." );
		return false;
	}
	
	event.SetInt( 50 );
	InitializeDTIColorMap( wnd, event );
	
	event.SetInt( 100 );
	wxPostEvent( wnd, event );
	
	return true;	
}

void LayerDTI::InitializeDTIColorMap( wxWindow* wnd, wxCommandEvent& event )
{
	vtkImageData* rasDTI = m_vectorSource->GetImageOutput();
	int* dim = rasDTI->GetDimensions();
	int nSize = dim[0]*dim[1]*dim[2];
	double v[3];
	int c[3];
	vtkDataArray* vectors = rasDTI->GetPointData()->GetScalars();
	vtkFloatArray* fas = vtkFloatArray::New();
	fas->DeepCopy( m_imageData->GetPointData()->GetScalars() );
	m_imageData->SetNumberOfScalarComponents( 2 );
	m_imageData->AllocateScalars();
	int nProgressStep = ( 99-event.GetInt() ) / 5;	
	for ( int i = 0; i < nSize; i++ )
	{
		vectors->GetTuple( i, v );
		double fa = fas->GetComponent( i, 0 );
		for ( int j = 0; j < 3; j++ )
		{
			c[j] = (int)(fabs(v[j]) * fa * 64);
			if ( c[j] > 63 )
				c[j] = 63;
		}
		float scalar = c[0]*64*64 + c[1]*64 + c[2];
		int x = i%dim[0];
		int y = (i/dim[0])%dim[1];
		int z = i/(dim[0]*dim[1]);
		m_imageData->SetScalarComponentFromFloat( x, y, z, 0, fa );
		m_imageData->SetScalarComponentFromFloat( x, y, z, 1, scalar );
		if ( nSize >= 5 && i%(nSize/5) == 0 )
		{
			event.SetInt( event.GetInt() + nProgressStep );  
			wxPostEvent( wnd, event );
		}
	}

	fas->Delete();
}

void LayerDTI::UpdateColorMap()
{
	if ( mProperties->GetColorMap() == LayerPropertiesMRI::DirectionCoded )
	{
		for ( int i = 0; i < 3; i++ )
		{
			mColorMap[i]->SetLookupTable( ( (LayerPropertiesDTI*)mProperties )->GetDirectionCodedTable() );
			mColorMap[i]->SetActiveComponent( 1 );
		}
	}
	else
		LayerMRI::UpdateColorMap();
}
