/**
 * @file  Contour2D.cpp
 * @brief Contour2D.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2010/03/09 19:47:23 $
 *    $Revision: 1.2 $
 *
 * Copyright (C) 2008-2009,
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

#include "Contour2D.h"
#include "RenderView2D.h"
#include "vtkImageData.h"
#include "vtkImageActor.h"
#include "vtkImageResample.h"
#include "vtkImageThreshold.h"
#include "vtkSimpleLabelEdgeFilter.h"
#include "vtkImageMapToColors.h"
#include "vtkRGBAColorTransferFunction.h"
#include "vtkMatrix4x4.h"
#include "vtkImageMask.h"

Contour2D::Contour2D( RenderView2D* view ) :
    Broadcaster( "Contour2D" ),
    Listener( "Contour2D" ),
    m_view( view )
{
  m_nPlane = view->GetViewPlane();
  m_actorContour = vtkSmartPointer<vtkImageActor>::New();
  m_actorContour->VisibilityOff();
  m_actorContour->InterpolateOff();
  
  m_filterThreshold = vtkSmartPointer<vtkImageThreshold>::New();
  m_filterThreshold->SetOutputScalarTypeToUnsignedChar();
  m_filterThreshold->ReplaceInOn();
  m_filterThreshold->ReplaceOutOn();
  m_filterThreshold->SetInValue( 1 );
  m_filterThreshold->SetOutValue( 0 );
  m_filterMask = vtkSmartPointer<vtkImageMask>::New();
  m_filterResample = vtkSmartPointer<vtkImageResample>::New();
  m_filterResample->SetAxisMagnificationFactor( 0, 2.0 );
  m_filterResample->SetAxisMagnificationFactor( 1, 2.0 );
  m_filterResample->SetAxisMagnificationFactor( 2, 2.0 );
  m_filterResample->SetInterpolationModeToNearestNeighbor();
  m_filterEdge = vtkSmartPointer<vtkSimpleLabelEdgeFilter>::New();
  m_colormap = vtkSmartPointer<vtkImageMapToColors>::New();
//  m_colormap->SetLookupTable( GetProperties()->GetGrayScaleTable() );
  m_colormap->SetOutputFormatToRGBA();
  m_colormap->PassAlphaToOutputOn();
  vtkSmartPointer<vtkRGBAColorTransferFunction> lut = vtkSmartPointer<vtkRGBAColorTransferFunction>::New();
  lut->AddRGBAPoint( 0, 0, 0, 0, 0 );
  lut->AddRGBAPoint( 1, 1, 1, 1, 1 );
  lut->Build();
  m_colormap->SetLookupTable( lut );
}

Contour2D::~Contour2D()
{}

vtkImageActor* Contour2D::GetActor()
{
  return m_actorContour.GetPointer();
}

void Contour2D::Reset()
{
}

vtkImageData* Contour2D::GetThresholdedImage()
{
  return m_filterMask->GetOutput();
}

void Contour2D::SetInput( vtkImageData* imagedata, double dContourValue, double dSliceLocation )
{
  m_filterThreshold->SetInput( imagedata );
  SetContourValue( dContourValue );
  m_filterMask->SetImageInput( m_filterThreshold->GetOutput() );
  
  // create a mask and initialize it.
  m_imageMask = vtkSmartPointer<vtkImageData>::New();
  m_imageMask->DeepCopy( m_filterThreshold->GetOutput() );
  int* dim = m_imageMask->GetDimensions();
  int size = dim[0]*dim[1]*dim[2];
  unsigned char* ptr = (unsigned char*)m_imageMask->GetScalarPointer();
  for ( int i = 0; i < size; i++ )
    ptr[i] = 1;
  
  m_filterMask->SetMaskInput( m_imageMask );
  m_filterResample->SetInput( m_filterMask->GetOutput() );
  m_filterEdge->SetInput( m_filterResample->GetOutput() );
  m_colormap->SetInput( m_filterEdge->GetOutput() );
  m_actorContour->SetInput( m_colormap->GetOutput() ); 
  
  UpdateSliceLocation( dSliceLocation );
}

void Contour2D::AddPatchLineOnMask( double* ras1, double* ras2 )
{
  if ( m_imageMask.GetPointer() == NULL )
    return;
  
  int n1[3], n2[3];
  double* origin = m_imageMask->GetOrigin();
  double* vsize = m_imageMask->GetSpacing();
  int* dim = m_imageMask->GetDimensions();
  for ( int i = 0; i < 3; i++ )
  {
    n1[i] = (int)( (ras1[i] - origin[i]) / vsize[i] );
    n2[i] = (int)( (ras2[i] - origin[i]) / vsize[i] );
  }
  
  int nx = 0, ny = 1;
  switch ( m_nPlane )
  {
    case 0:
      nx = 1;
      ny = 2;
      break;
    case 1:
      nx = 0;
      ny = 2;
      break;
  }
  
  unsigned char* ptr = (unsigned char*)m_imageMask->GetScalarPointer();
  int x0 = n1[nx], y0 = n1[ny], x1 = n2[nx], y1 = n2[ny];
  int dx = x1 - x0;
  int dy = y1 - y0;
  double t = 0.5;
  int n[3];
  if ( m_nPlane == 0 )
    ptr[n1[ny]*dim[0]+dim[0]-n1[nx]-1] = 0;
  else
    ptr[n1[ny]*dim[0]+n1[nx]] = 0;
  if ( abs( dx ) > abs( dy ) )
  {
    double m = (double) dy / (double) dx;
    t += y0;
    dx = ( dx < 0 ? -1 : 1 );
    m *= dx;
    while ( x0 != x1 )
    {
      x0 += dx;
      t += m;
      n[nx] = x0;
      n[ny] = (int) t;
      if ( m_nPlane == 0 )
        ptr[n[ny]*dim[0]+dim[0]-n[nx]-1] = 0;
      else
        ptr[n[ny]*dim[0]+n[nx]] = 0;
    }
  }
  else
  {
    double m = (double) dx / (double) dy;
    t += x0;
    dy = ( dy < 0 ? -1 : 1 );
    m *= dy;
    while ( y0 != y1 )
    {
      y0 += dy;
      t += m;
      n[nx] = (int) t;
      n[ny] = y0;
      if ( m_nPlane == 0 )
        ptr[n[ny]*dim[0]+dim[0]-n[nx]-1] = 0;
      else
        ptr[n[ny]*dim[0]+n[nx]] = 0;
    }
  }
  
  m_imageMask->Modified();
}

void Contour2D::UpdateSliceLocation( double dSliceLocation )
{
  vtkImageData* imagedata = vtkImageData::SafeDownCast( m_filterThreshold->GetInput() );
  if ( !imagedata )
    return;
  
  if ( fabs( dSliceLocation - m_dSliceLocation ) < 1e-6 )
    return;
  
  m_dSliceLocation = dSliceLocation;
  vtkSmartPointer<vtkMatrix4x4> matrix =
      vtkSmartPointer<vtkMatrix4x4>::New();
  matrix->Identity();
  double* vsize = imagedata->GetSpacing();
  double pos[3] = { vsize[0]/4.0, vsize[1]/4.0, vsize[2]/4.0 };
  switch ( m_nPlane )
  {
    case 0:
      m_actorContour->PokeMatrix( matrix );
      m_actorContour->SetPosition( dSliceLocation, -pos[1], pos[2] );
      m_actorContour->RotateX( 90 );
      m_actorContour->RotateY( -90 );
      break;
    case 1:
      m_actorContour->PokeMatrix( matrix );
      m_actorContour->SetPosition( pos[0], dSliceLocation, pos[2] );
      m_actorContour->RotateX( 90 );
      break;
    case 2:
      m_actorContour->SetPosition( pos[0], pos[1], dSliceLocation );
      break;
  }
}

void Contour2D::SetContourValue( double dContourValue )
{
  m_dContourValue = dContourValue;
  m_filterThreshold->ThresholdByUpper( dContourValue );
  m_filterThreshold->Update();
}

bool Contour2D::IsVisible()
{
  return m_actorContour->GetVisibility();
}

void Contour2D::SetVisible( bool visible )
{
  m_actorContour->SetVisibility( visible?1:0 );
}
