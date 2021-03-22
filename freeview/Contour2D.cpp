/**
 * @brief Contour2D.
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

#include "Contour2D.h"
#include "RenderView2D.h"
#include "vtkImageData.h"
#include "vtkImageActor.h"
#include "vtkImageReslice.h"
#include "vtkImageThreshold.h"
#include "vtkSimpleLabelEdgeFilter.h"
#include "vtkImageMapToColors.h"
#include "vtkRGBAColorTransferFunction.h"
#include "vtkImageGaussianSmooth.h"
#include "vtkImageExtractComponents.h"
#include "vtkMatrix4x4.h"
#include "vtkImageMask.h"
#include "vtkImageLogic.h"
#include "vtkImageMapper3D.h"
#include <QDebug>

#define IMAGE_RESAMPLE_FACTOR     4.0     // must be multiples of 2

Contour2D::Contour2D( RenderView2D* view ) :
  QObject( view ),
  m_view( view )
{
  m_nPlane = view->GetViewPlane();
  m_imageInput = vtkSmartPointer<vtkImageData>::New();
  m_dContourValue = 0;

  m_actorContour = vtkSmartPointer<vtkImageActor>::New();
  m_actorContour->VisibilityOff();
  m_actorContour->InterpolateOff();
#if VTK_MAJOR_VERSION > 5
    m_actorContour->ForceOpaqueOn();
#endif

  m_filterSmooth = vtkSmartPointer<vtkImageGaussianSmooth>::New();
  m_filterSmooth->SetStandardDeviations( 1, 1, 1 );
  m_filterThreshold = vtkSmartPointer<vtkImageThreshold>::New();
  m_filterThreshold->SetOutputScalarTypeToUnsignedChar();
  m_filterThreshold->ReplaceInOn();
  m_filterThreshold->ReplaceOutOn();
  m_filterThreshold->SetInValue( 1 );
  m_filterThreshold->SetOutValue( 0 );
  m_filterMask  = vtkSmartPointer<vtkImageMask>::New();
  m_filterLogic = vtkSmartPointer<vtkImageLogic>::New();
  m_filterLogic->SetOperationToOr();
  m_filterResample = vtkSmartPointer<vtkImageReslice>::New();
//  m_filterResample->SetAxisMagnificationFactor( 0, IMAGE_RESAMPLE_FACTOR );
//  m_filterResample->SetAxisMagnificationFactor( 1, IMAGE_RESAMPLE_FACTOR );
//  m_filterResample->SetAxisMagnificationFactor( 2, IMAGE_RESAMPLE_FACTOR );
  m_filterResample->SetInterpolationModeToNearestNeighbor();
  m_filterEdge = vtkSmartPointer<vtkSimpleLabelEdgeFilter>::New();
  m_colormap = vtkSmartPointer<vtkImageMapToColors>::New();
  m_colormap->SetOutputFormatToRGBA();
  m_colormap->PassAlphaToOutputOn();

  SetColor( 1, 1, 1 );

  connect(this, SIGNAL(ValueChanged()), view, SLOT(RequestRedraw()));
  connect(this, SIGNAL(ColorChanged()), view, SLOT(RequestRedraw()));

  Reset();
}

Contour2D::~Contour2D()
{}

vtkImageActor* Contour2D::GetActor()
{
  return m_actorContour;
}

vtkImageData* Contour2D::GetInputImage()
{
  return m_imageInput;
}

void Contour2D::Reset()
{
  m_imageInput = NULL;
}

vtkImageData* Contour2D::GetThresholdedImage()
{
  if (m_filterMask->GetInput())
    return m_filterMask->GetOutput();
  else
    return NULL;
}

void Contour2D::SetInput( vtkImageData* imagedata, double dContourValue, double dSliceLocation, int active_frame )
{
  vtkSmartPointer<vtkImageExtractComponents> extract = vtkSmartPointer<vtkImageExtractComponents>::New();
  if ( imagedata->GetNumberOfScalarComponents() > 1 )
  {
#if VTK_MAJOR_VERSION > 5
    extract->SetInputData( imagedata );
#else
    extract->SetInput(imagedata);
#endif
    extract->SetComponents( active_frame );
    extract->Update();
    m_imageInput = extract->GetOutput();
  }
  else
  {
    m_imageInput = imagedata;
  } 

#if VTK_MAJOR_VERSION > 5
  m_filterThreshold->SetInputData( m_imageInput );
#else
  m_filterThreshold->SetInput(m_imageInput);
#endif
  SetContourValue( dContourValue );

  // create two masks and initialize them.
  m_imageMaskAdd = vtkSmartPointer<vtkImageData>::New();
  m_imageMaskRemove = vtkSmartPointer<vtkImageData>::New();
  m_imageMaskAdd->DeepCopy( m_filterThreshold->GetOutput() );
  m_imageMaskRemove->DeepCopy( m_imageMaskAdd );
  int* dim = m_imageMaskAdd->GetDimensions();
  long long size = ((long long)dim[0])*dim[1]*dim[2];
  memset( m_imageMaskAdd->GetScalarPointer(), 0, size );
  unsigned char* ptr = (unsigned char*)m_imageMaskRemove->GetScalarPointer();
  for ( int i = 0; i < size; i++ )
  {
    ptr[i] = 1;
  }

#if VTK_MAJOR_VERSION > 5
  m_filterLogic->SetInput1Data( m_filterThreshold->GetOutput() );
  m_filterLogic->SetInput2Data( m_imageMaskAdd );
  m_filterMask->SetMaskInputData( m_imageMaskRemove );
#else
  m_filterLogic->SetInput1( m_filterThreshold->GetOutput() );
  m_filterLogic->SetInput2( m_imageMaskAdd );
  m_filterMask->SetMaskInput( m_imageMaskRemove );
#endif
  m_filterMask->SetInputConnection( m_filterLogic->GetOutputPort() );
  double vs[3];
  imagedata->GetSpacing(vs);
  m_filterResample->SetOutputSpacing(vs[0]/IMAGE_RESAMPLE_FACTOR, vs[1]/IMAGE_RESAMPLE_FACTOR, vs[2]/IMAGE_RESAMPLE_FACTOR);
  m_filterResample->SetInputConnection( m_filterMask->GetOutputPort() );
  m_filterEdge->SetInputConnection( m_filterResample->GetOutputPort() );
  m_colormap->SetInputConnection( m_filterEdge->GetOutputPort() );
#if VTK_MAJOR_VERSION > 5
  m_actorContour->GetMapper()->SetInputConnection( m_colormap->GetOutputPort() );
#else
  m_actorContour->SetInput( m_colormap->GetOutput() );
#endif

  SetSmooth( m_bSmooth );
  UpdateSliceLocation( dSliceLocation, true );
}

void Contour2D::AddLine( double* ras1, double* ras2 )
{
  DrawPatchLineOnMask( m_imageMaskAdd, ras1, ras2, 1 );
  // update the remove mask as well
  DrawPatchLineOnMask( m_imageMaskRemove, ras1, ras2, 1 );
}

void Contour2D::RemoveLine( double* ras1, double* ras2 )
{
  DrawPatchLineOnMask( m_imageMaskRemove, ras1, ras2, 0 );
}

void Contour2D::DrawPatchLineOnMask( vtkImageData* image, double* ras1, double* ras2, int nDrawValue )
{
  if ( !image )
  {
    return;
  }

  int n1[2], n2[2];
  double* origin = image->GetOrigin();    // 2D image!
  double* vsize = image->GetSpacing();
  int* dim = image->GetDimensions();

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
  n1[0] = (int)( (ras1[nx] - origin[0]) / vsize[0] + 0.5 );
  n1[1] = (int)( (ras1[ny] - origin[1]) / vsize[1] + 0.5 );
  n2[0] = (int)( (ras2[nx] - origin[0]) / vsize[0] + 0.5 );
  n2[1] = (int)( (ras2[ny] - origin[1]) / vsize[1] + 0.5 );

  nx = 0;
  ny = 1;
  unsigned char* ptr = (unsigned char*)image->GetScalarPointer();
  int x0 = n1[nx], y0 = n1[ny], x1 = n2[nx], y1 = n2[ny];
  int dx = x1 - x0;
  int dy = y1 - y0;
  double t = 0.5;
  int n[2];
  ptr[n1[ny]*dim[0]+n1[nx]] = nDrawValue;
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
      ptr[n[ny]*dim[0]+n[nx]] = nDrawValue;
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
      ptr[n[ny]*dim[0]+n[nx]] = nDrawValue;
    }
  }

  image->Modified();
}

void Contour2D::UpdateSliceLocation( double dSliceLocation, bool bForced )
{
  vtkImageData* imagedata = vtkImageData::SafeDownCast( m_filterThreshold->GetInput() );
  if ( !imagedata )
  {
    return;
  }

  if ( !bForced && fabs( dSliceLocation - m_dSliceLocation ) < 1e-6 )
  {
    return;
  }

  m_dSliceLocation = dSliceLocation;
  vtkSmartPointer<vtkMatrix4x4> matrix =
      vtkSmartPointer<vtkMatrix4x4>::New();
  matrix->Identity();
  double* vsize = imagedata->GetSpacing();    // 2D spacing!
  double pos[2] = { vsize[0]/IMAGE_RESAMPLE_FACTOR/2, vsize[1]/IMAGE_RESAMPLE_FACTOR/2 };
  switch ( m_nPlane )
  {
  case 0:
    m_actorContour->PokeMatrix( matrix );
    m_actorContour->SetPosition( dSliceLocation, pos[0], pos[1] );
    m_actorContour->RotateX( 90 );
    m_actorContour->RotateY( 90 );
    break;
  case 1:
    m_actorContour->PokeMatrix( matrix );
    m_actorContour->SetPosition( pos[0], dSliceLocation, pos[1] );
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
  emit ValueChanged();
}

bool Contour2D::IsVisible()
{
  return m_actorContour->GetVisibility();
}

void Contour2D::SetVisible( bool visible )
{
  m_actorContour->SetVisibility( visible?1:0 );
}

void Contour2D::SetSmooth( bool bSmooth )
{
  m_bSmooth = bSmooth;
  if ( m_imageInput.GetPointer() )
  {
#if VTK_MAJOR_VERSION > 5
    m_filterSmooth->SetInputData( m_imageInput );
    m_filterThreshold->SetInputData( bSmooth ? m_filterSmooth->GetOutput() : m_imageInput.GetPointer() );
#else
    m_filterSmooth->SetInput( m_imageInput );
    m_filterThreshold->SetInput( bSmooth ? m_filterSmooth->GetOutput() : m_imageInput.GetPointer() );
#endif
  }
}

double Contour2D::GetSmoothSD()
{
  double* sd = m_filterSmooth->GetStandardDeviations();
  return sd[0];
}

void Contour2D::SetSmoothSD( double sd )
{
  m_filterSmooth->SetStandardDeviations( sd, sd, sd );
}

void Contour2D::SetColor( double r, double g, double b )
{
  m_dContourColor[0] = r;
  m_dContourColor[1] = g;
  m_dContourColor[2] = b;

  vtkSmartPointer<vtkRGBAColorTransferFunction> lut = vtkSmartPointer<vtkRGBAColorTransferFunction>::New();
  lut->AddRGBAPoint( 0, 0, 0, 0, 0 );
  lut->AddRGBAPoint( 1, r, g, b, 1 );
  lut->Build();
  m_colormap->SetLookupTable( lut );

  if ( IsVisible() )
  {
    emit ColorChanged();
  }
}

