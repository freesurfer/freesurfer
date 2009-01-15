/**
 * @file  LivewireTool.cpp
 * @brief LivewireTool.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2009/01/15 19:28:58 $
 *    $Revision: 1.3 $
 *
 * Copyright (C) 2002-2009,
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
 
#include "LivewireTool.h"
#include <vtkPoints.h>
#include <vtkImageClip.h>
#include <vtkImageData.h>
#include <vtkImageGradientMagnitude.h>
#include <vtkImageShiftScale.h>
#include <vtkImageChangeInformation.h>
#include <vtkImageAnisotropicDiffusion2D.h>
#include <vtkImageGaussianSmooth.h>
#include <vtkCellArray.h>
#include "vtkDijkstraImageGeodesicPath.h"

LivewireTool::LivewireTool( ) :
	m_nPlane( 0 ),
	m_nSlice( 0 ),
	m_imageData( NULL ),
	m_imageSlice( NULL )
{
	m_imageClip = vtkSmartPointer<vtkImageClip>::New();
	m_path = vtkSmartPointer<vtkDijkstraImageGeodesicPath>::New();
	m_info = vtkSmartPointer<vtkImageChangeInformation>::New();
}

LivewireTool::~LivewireTool()
{
}

void LivewireTool::SetImagePlane( int nPlane )
{ 
	if ( m_imageData && m_nPlane != nPlane )
	{
		m_nPlane = nPlane; 	
		int ext[6];
		m_imageData->GetExtent( ext );
		ext[m_nPlane*2] = ext[m_nPlane*2 + 1] = m_nSlice;
		m_imageClip->SetOutputWholeExtent( ext );
		int n[3] = { 0, 0, 0 };
		n[m_nPlane] = -1*m_nSlice;
		m_info->SetExtentTranslation( n );
	//	m_info->SetOutputOrigin( 0, 0, 0 );
		m_path->Update();
	}
}
			
void LivewireTool::SetImageSlice( int nSlice )
{ 
	if ( m_imageData && m_nSlice != nSlice )
	{
		m_nSlice = nSlice; 
		int ext[6];
		m_imageData->GetExtent( ext );
		ext[m_nPlane*2] = ext[m_nPlane*2 + 1] = m_nSlice;
		m_imageClip->SetOutputWholeExtent( ext );
		int n[3] = { 0, 0, 0 };
		n[m_nPlane] = -1*m_nSlice;
		m_info->SetExtentTranslation( n );
	//	m_info->SetOutputOrigin( 0, 0, 0 );
		m_path->Update();
	//	m_imageSlice->GetExtent( ext );
	//	cout << ext[0] << " " << ext[1] <<  " " << ext[2] << " " << ext[3] << " " << ext[4] << " " << ext[5] << endl;
		double* orig = m_info->GetOutput()->GetOrigin();
	//	cout << "orig " << orig[0] << " " << orig[1] << " " << orig[2] << endl;
		orig = m_imageData->GetOrigin();
	//	cout << "orig " << orig[0] << " " << orig[1] << " " << orig[2] << endl;
		m_imageSlice->GetExtent( ext );
	//	cout << "ext " << ext[0] << " " << ext[1] << " " << ext[2] << " " << ext[3] << " " << ext[4] << " " << ext[5] << endl;	
	}
}

void LivewireTool::SetImageData( vtkImageData* image_in )
{
	if ( m_imageData != image_in )
	{
		m_imageData = image_in;
		vtkSmartPointer<vtkImageChangeInformation> info = vtkSmartPointer<vtkImageChangeInformation>::New();
		info->SetOutputOrigin( 0, 0, 0 );
		info->SetOutputExtentStart( 0, 0, 0 );
		info->SetInput( image_in );
		info->Update();
		double* orig = info->GetOutput()->GetOrigin();
		int ext[6];
		info->GetOutput()->GetExtent(ext);
		cout << "orig " << orig[0] << " " << orig[1] << " " << orig[2] << endl;
		cout << "ext " << ext[0] << " " << ext[1] << " " << ext[2] << " " << ext[3] << " " << ext[4] << " " << ext[5] << endl;	
	
		m_imageClip->SetInput( info->GetOutput() );
	//	int ext[6];
		image_in->GetExtent( ext );
		ext[m_nPlane*2] = ext[m_nPlane*2 + 1] = m_nSlice;
		m_imageClip->SetOutputWholeExtent( ext );
//		m_imageClip->ClipDataOn();
		m_imageClip->ReleaseDataFlagOff();
		m_imageClip->Update();
		
		vtkSmartPointer<vtkImageAnisotropicDiffusion2D> smooth = vtkSmartPointer<vtkImageAnisotropicDiffusion2D>::New();
		smooth->SetInputConnection( m_imageClip->GetOutputPort() );
		smooth->SetDiffusionFactor( 0.75 );
		smooth->SetDiffusionThreshold( 50.0 );
		smooth->SetNumberOfIterations( 5 ); 
	
/*	vtkSmartPointer<vtkImageGaussianSmooth> smooth = vtkSmartPointer<vtkImageGaussianSmooth>::New();
		smooth->SetInputConnection( clip->GetOutputPort() );
		smooth->SetStandardDeviations( 1, 1, 1 );*/
	
		vtkSmartPointer<vtkImageGradientMagnitude> grad = vtkSmartPointer<vtkImageGradientMagnitude>::New();
		grad->SetDimensionality( 2 );
		grad->HandleBoundariesOn();
		grad->SetInputConnection( smooth->GetOutputPort() );
		grad->Update();		
		
		double* range = grad->GetOutput()->GetScalarRange();
		vtkSmartPointer<vtkImageShiftScale> scale = vtkSmartPointer<vtkImageShiftScale>::New();
		scale->SetShift( -1.0*range[1] );
		scale->SetScale( 255.0 /( range[0] - range[1] ) );
		scale->SetOutputScalarTypeToShort();
		scale->SetInputConnection( grad->GetOutputPort() );
		scale->ReleaseDataFlagOff();

		m_info->SetInput( scale->GetOutput() );
		int n[3] = { 0, 0, 0 };
		n[m_nPlane] = -1*m_nSlice;
		m_info->SetExtentTranslation( n );
		m_info->Update();

		m_imageSlice = scale->GetOutput();
		m_path->SetInput( m_info->GetOutput() );
	//	m_path->Update();
	}
}
		
void LivewireTool::GetLivewirePoints( double* pt1_in, double* pt2_in, vtkPoints* pts_out )
{
	if ( !m_imageSlice )
		return; 
	
	double pt1[3], pt2[3];
	double* orig = m_imageData->GetOrigin();
	for ( int i = 0; i < 3; i++ )
	{
		pt1[i] = pt1_in[i] - orig[i];
		pt2[i] = pt2_in[i] - orig[i];
	}
	
	vtkIdType beginVertId = m_imageSlice->FindPoint( pt1 );
	vtkIdType endVertId = m_imageSlice->FindPoint( pt2 );
	cout << beginVertId << "  " << endVertId << endl;

	if ( beginVertId == -1 || endVertId == -1 ) 
	{
	//	cout << "can not find point: " << pt1_in[0] << " " << pt1_in[1] << " " << pt1_in[2] << ", " 
	//			<< pt2_in[0] << " " << pt2_in[1] << " " << pt2_in[2] << endl;
		return;
	}

	m_path->SetStartVertex( endVertId );
	m_path->SetEndVertex( beginVertId );
	m_path->Update();
	double* pt = m_info->GetOutput()->GetPoint( beginVertId );
	cout << pt[0] << " " << pt[1] << " " << pt[2] << endl;

	vtkPolyData *pd = m_path->GetOutput();
	vtkIdType npts = 0, *pts = NULL;
	pd->GetLines()->InitTraversal();
	pd->GetLines()->GetNextCell( npts, pts );
	cout << npts << endl;
	double offset[3] = { 0, 0, 0 };
	double* vs = m_imageData->GetSpacing();
	offset[m_nPlane] = m_nSlice*vs[m_nPlane];
	for ( int i = 0; i < npts; i++ )
	{
		double* p = pd->GetPoint( pts[i] );
	//	cout << p[0] << " " << p[1] << " " << p[2] << endl;
		pts_out->InsertNextPoint( p[0] + offset[0], p[1] + offset[1], p[2] + offset[2] );
	}
}


