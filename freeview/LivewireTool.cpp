/**
 * @brief LivewireTool.
 *
 */
/*
 * Original Author: Ruopeng Wang
 *
 * Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
 *
 * Terms and conditions for use, reproduction, distribution and contribution
 * are found in the 'FreeSurfer Software License Agreement' contained
 * in the file 'LICENSE' found in the FreeSurfer distribution, and here:
 *
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
 *
 * Reporting: freesurfer@nmr.mgh.harvard.edu
 *
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
  m_path = vtkSmartPointer<vtkDijkstraImageGeodesicPath>::New();
}

LivewireTool::~LivewireTool()
{}

void LivewireTool::GetLivewirePoints( vtkImageData* image, int nPlane, int nSlice,
                                      double* pt1_in, double* pt2_in, vtkPoints* pts_out )
{
  UpdateImageDataInfo( image, nPlane, nSlice );
  GetLivewirePoints( pt1_in, pt2_in, pts_out );
}

void LivewireTool::UpdateImageDataInfo( vtkImageData* image_in, int nPlane, int nSlice )
{
  if ( m_imageData != image_in || m_nPlane != nPlane || m_nSlice != nSlice )
  {
    m_imageData = image_in;
    m_nPlane = nPlane;
    m_nSlice = nSlice;
    vtkSmartPointer<vtkImageClip> clip = vtkSmartPointer<vtkImageClip>::New();
#if VTK_MAJOR_VERSION > 5
    clip->SetInputData( image_in );
#else
    clip->SetInput( image_in );
#endif
    int ext[6];
    image_in->GetExtent( ext );
    ext[m_nPlane*2] = ext[m_nPlane*2 + 1] = m_nSlice;
    clip->SetOutputWholeExtent( ext );
    clip->ClipDataOn();
    clip->ReleaseDataFlagOff();

    vtkSmartPointer<vtkImageAnisotropicDiffusion2D> smooth = vtkSmartPointer<vtkImageAnisotropicDiffusion2D>::New();
    smooth->SetInputConnection( clip->GetOutputPort() );
    smooth->SetDiffusionFactor( 0.75 );
    smooth->SetDiffusionThreshold( 50.0 );
    smooth->SetNumberOfIterations( 5 );

    /*  vtkSmartPointer<vtkImageGaussianSmooth> smooth = vtkSmartPointer<vtkImageGaussianSmooth>::New();
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

    vtkSmartPointer<vtkImageChangeInformation> info = vtkSmartPointer<vtkImageChangeInformation>::New();
    info->SetInputConnection( scale->GetOutputPort() );
    int n[3] = { 0, 0, 0 };
    n[m_nPlane] = -1*m_nSlice;
    info->SetExtentTranslation( n );
    info->Update();

    m_imageSlice = scale->GetOutput();
    m_path->SetInputConnection( info->GetOutputPort() );
  }
}

void LivewireTool::GetLivewirePoints( double* pt1_in, double* pt2_in, vtkPoints* pts_out )
{
  if ( !m_imageSlice )
  {
    return;
  }

  vtkIdType beginVertId = m_imageSlice->FindPoint( pt1_in );
  vtkIdType endVertId = m_imageSlice->FindPoint( pt2_in );
  // cout << beginVertId << "  " << endVertId << endl;

  if ( beginVertId == -1 || endVertId == -1 )
  {
    // cout << "can not find point: " << pt1_in[0] << " " << pt1_in[1] << " " << pt1_in[2] << ", "
    //   << pt2_in[0] << " " << pt2_in[1] << " " << pt2_in[2] << endl;
    return;
  }

  m_path->SetStartVertex( endVertId );
  m_path->SetEndVertex( beginVertId );
  m_path->Update();
  // double* pt = m_info->GetOutput()->GetPoint( beginVertId );
  // cout << pt[0] << " " << pt[1] << " " << pt[2] << endl;

  vtkPolyData *pd = m_path->GetOutput();
  vtkIdType npts = 0, *pts = NULL;
  pd->GetLines()->InitTraversal();
  pd->GetLines()->GetNextCell( npts, pts );
  // cout << npts << endl;
  double offset[3] = { 0, 0, 0 };
  double* vs = m_imageData->GetSpacing();
  offset[m_nPlane] = m_nSlice*vs[m_nPlane];
  for ( int i = 0; i < npts; i++ )
  {
    double* p = pd->GetPoint( pts[i] );
    // cout << p[0] << " " << p[1] << " " << p[2] << endl;
    pts_out->InsertNextPoint( p[0] + offset[0], p[1] + offset[1], p[2] + offset[2] );
  }
}


