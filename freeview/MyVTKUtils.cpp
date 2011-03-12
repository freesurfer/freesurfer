/**
 * @file  MyVTKUtilsUtils.h
 * @brief Misc utility class.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: krish $
 *    $Date: 2011/03/12 00:28:52 $
 *    $Revision: 1.1 $
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

#include "MyVTKUtils.h"
#include <math.h>
#include <stddef.h>
#include <vtkSmartPointer.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkImageWriter.h>
#include <vtkBMPWriter.h>
#include <vtkJPEGWriter.h>
#include <vtkPNGWriter.h>
#include <vtkTIFFWriter.h>
#include <vtkPostScriptWriter.h>
#include <vtkRenderLargeImage.h>
#include <vtkVRMLExporter.h>
#include <vtkActor.h>
#include <vtkImageData.h>
#include <vtkImageDilateErode3D.h>
#include <vtkContourFilter.h>
#include <vtkMarchingContourFilter.h>
#include <vtkPolyDataConnectivityFilter.h>
#include <vtkPolyData.h>
#include <vtkCellArray.h>
#include <vtkSmoothPolyDataFilter.h>
#include <vtkPolyDataMapper.h>
#include <vtkStripper.h>
#include <vtkPolyDataNormals.h>
#include <vtkCamera.h>
#include <vtkCubeSource.h>
#include <vtkLineSource.h>
#include <vtkSphereSource.h>
#include <vtkTubeFilter.h>
#include <vtkProperty.h>
#include <vtkVectorText.h>
#include <vtkFollower.h>
#include <vtkAppendPolyData.h>
#include <vtkConeSource.h>
#include <vtkImageThreshold.h>
#include <vtkDecimatePro.h>
#include <vtkTriangleFilter.h>
#include <vtkPiecewiseFunction.h>
#include <vtkColorTransferFunction.h>
#include <vtkVolumeRayCastMapper.h>
#include <vtkVolumeProperty.h>
#include <vtkVolume.h>
#include <vtkVolumeRayCastCompositeFunction.h>
#include <vtkVolumeTextureMapper2D.h>
#include <vtkFixedPointVolumeRayCastMapper.h>
#include <vtkImageCast.h>
#include <vtkImageClip.h>
#include <vtkImageGradientMagnitude.h>
#include <vtkImageShiftScale.h>
#include <vtkImageChangeInformation.h>
#include <vtkImageAnisotropicDiffusion2D.h>
#include <vtkImageGaussianSmooth.h>
#include <vtkImageClip.h>
#include <vtkDijkstraImageGeodesicPath.h>
#include <QFileInfo>

bool MyVTKUtils::VTKScreenCapture( vtkRenderWindow* renderWnd,
                                vtkRenderer* renderer,
                                const char* filename,
                                bool bAntiAliasing,
                                int nMag )
{
  QString fn = filename;
  vtkImageWriter* writer = 0;
  QString ext = QFileInfo(filename).suffix();
  if ( ext == "wrl" )
  {
    vtkVRMLExporter* exporter = vtkVRMLExporter::New();
    exporter->SetFileName( filename );
    exporter->SetRenderWindow( renderWnd );
    exporter->Write();
    exporter->Delete();
  }
  else if ( ext == "jpg" || ext == "jpeg" )
    writer = vtkJPEGWriter::New();
  else if ( ext == "bmp" )
    writer = vtkBMPWriter::New();
  else if ( ext == "ps" )
    writer = vtkPostScriptWriter::New();
  else if ( ext == "tif" || ext == "tiff" )
    writer = vtkTIFFWriter::New();
  else
  {
    writer = vtkPNGWriter::New();
    if ( ext != "png" )
      fn += ".png";
  }

  bool ret = true;
  if (writer)
  {
    // bool bCurrentAA = GetAntialiasing() > 0;
    // SetAntialiasing(bAntiAliasing, false);
    vtkRenderLargeImage* image = vtkRenderLargeImage::New();
    image->SetInput( renderer );
    image->SetMagnification( nMag );
    writer->SetInput( image->GetOutput() );
    writer->SetFileName( fn.toUtf8().data() );
    writer->Write();
    if ( writer->GetErrorCode() != 0 )
      ret = false;
    image->Delete();
    writer->Delete();
    // SetAntialiasing(bCurrentAA, false);
  }
  return ret;
}


void MyVTKUtils::ViewportToWorld( vtkRenderer* renderer, double x, double y, double z,
                               double& world_x, double& world_y, double& world_z )
{
  world_x = x;
  world_y = y;
  renderer->ViewportToNormalizedViewport( world_x, world_y );
  NormalizedViewportToWorld( renderer, world_x, world_y, z,
                             world_x, world_y, world_z );
}

void MyVTKUtils::ViewportToWorld( vtkRenderer* renderer,
                               double x, double y,
                               double& world_x,
                               double& world_y,
                               double& world_z )
{
  world_x = x;
  world_y = y;
  renderer->ViewportToNormalizedViewport( world_x, world_y );
  NormalizedViewportToWorld( renderer, world_x, world_y,
                             world_x, world_y, world_z );
}

void MyVTKUtils::NormalizedViewportToWorld( vtkRenderer* renderer,
                                         double x, double y, double z,
                                         double& world_x,
                                         double& world_y,
                                         double& world_z )
{
  world_x = x;
  world_y = y;
  world_z = z;
  renderer->NormalizedViewportToView( world_x, world_y, world_z );
  renderer->ViewToWorld( world_x, world_y, world_z );
}

void MyVTKUtils::NormalizedViewportToWorld( vtkRenderer* renderer,
                                         double x, double y,
                                         double& world_x,
                                         double& world_y,
                                         double& world_z )
{
  NormalizedViewportToWorld( renderer, x, y, 0.0, world_x, world_y, world_z );
}

void MyVTKUtils::WorldToViewport( vtkRenderer* renderer,
                               double world_x, double world_y, double world_z,
                               double& x, double& y, double& z )
{
  x = world_x;
  y = world_y;
  z = world_z;
  renderer->WorldToView( x, y, z );
  renderer->ViewToNormalizedViewport( x, y, z );
  renderer->NormalizedViewportToViewport( x, y );
}

/*
bool MyVTKUtils::BuildContourActor( vtkImageData* data_in,
                                 double dTh1, double dTh2,
                                 vtkActor* actor_out )
{
  vtkImageData* imagedata = data_in;

// int nValue = nThreshold;
  int nSwell = 2;
  vtkSmartPointer<vtkImageDilateErode3D> dilate =
    vtkSmartPointer<vtkImageDilateErode3D>::New();
  dilate->SetInput( imagedata );
  dilate->SetKernelSize( nSwell, nSwell, nSwell );
  dilate->SetDilateValue( dTh1 );
  dilate->SetErodeValue( 0 );
  vtkSmartPointer<vtkImageDilateErode3D> erode =
    vtkSmartPointer<vtkImageDilateErode3D>::New();
  erode->SetInput( dilate->GetOutput() );
  erode->SetKernelSize( 1, 1, 1 );
  erode->SetDilateValue( 0 );
  erode->SetErodeValue( dTh1 );

  vtkSmartPointer<vtkImageThreshold> threshold =
    vtkSmartPointer<vtkImageThreshold>::New();
  threshold->SetOutputScalarTypeToShort();
  threshold->SetInput( dilate->GetOutput() );
  threshold->ThresholdBetween( dTh1, dTh2+0.0001 );
  threshold->ReplaceOutOn();
  threshold->SetOutValue( 0 );
  vtkSmartPointer<vtkContourFilter> contour =
    vtkSmartPointer<vtkContourFilter>::New();
  contour->SetInput( threshold->GetOutput() );
  contour->SetValue( 0, dTh1 );
  vtkSmartPointer<vtkTriangleFilter> tri =
    vtkSmartPointer<vtkTriangleFilter>::New();
  tri->SetInput( contour->GetOutput() );
  vtkSmartPointer<vtkDecimatePro> decimate =
    vtkSmartPointer<vtkDecimatePro>::New();
  decimate->SetTargetReduction( 0.9 );
  decimate->SetInput( tri->GetOutput() );

  vtkPolyData* polydata = contour->GetOutput();
  polydata->Update();

  if ( polydata->GetNumberOfPoints() <= 0 )
  {
    return false;
  }
  else
  {
    vtkSmartPointer<vtkSmoothPolyDataFilter> smoother =
      vtkSmartPointer<vtkSmoothPolyDataFilter>::New();
    smoother->SetInput( contour->GetOutput() );
    smoother->SetNumberOfIterations( 30 );
    vtkSmartPointer<vtkPolyDataNormals> normals =
      vtkSmartPointer<vtkPolyDataNormals>::New();
    normals->SetInput( smoother->GetOutput()) ;
    normals->SetFeatureAngle( 90.0 );
    vtkSmartPointer<vtkStripper> stripper =
      vtkSmartPointer<vtkStripper>::New();
    stripper->SetInput( normals->GetOutput() );
    vtkSmartPointer<vtkPolyDataMapper> mapper =
      vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInput( decimate->GetOutput() );
    // mapper->ScalarVisibilityOff();
    actor_out->SetMapper( mapper );
    return true;
  }
}
*/

bool MyVTKUtils::BuildContourActor( vtkImageData* data_in,
                                 double dTh1, double dTh2,
                                 vtkActor* actor_out, int nSmoothIterations, int* ext, bool bAllRegions )
{
  double nValue = 1;
  int nSwell = 2;
  vtkSmartPointer<vtkImageThreshold> threshold = vtkSmartPointer<vtkImageThreshold>::New();

  if ( ext )
  {
    vtkSmartPointer<vtkImageClip> clipper = vtkSmartPointer<vtkImageClip>::New();
    clipper->SetInput( data_in );
    clipper->SetOutputWholeExtent( ext );
    threshold->SetInputConnection( clipper->GetOutputPort() );
  }
  else
    threshold->SetInput( data_in );
  threshold->ThresholdByLower( dTh2 );
  threshold->ReplaceOutOn();
  threshold->SetOutValue( dTh1-0.00001 );

  // dilate/erode is not used for now
  vtkSmartPointer<vtkImageDilateErode3D> dilate = vtkSmartPointer<vtkImageDilateErode3D>::New();
  dilate->SetInputConnection(threshold->GetOutputPort());
  dilate->SetKernelSize(nSwell, nSwell, nSwell);
  dilate->SetDilateValue(nValue);
  dilate->SetErodeValue(0);
  vtkSmartPointer<vtkImageDilateErode3D> erode = vtkSmartPointer<vtkImageDilateErode3D>::New();
  erode->SetInputConnection(dilate->GetOutputPort());
  erode->SetKernelSize(1, 1, 1);
  erode->SetDilateValue(0);
  erode->SetErodeValue(nValue);
  // end of dilate/erode

  vtkSmartPointer<vtkContourFilter> contour = vtkSmartPointer<vtkContourFilter>::New();
  contour->SetInputConnection( threshold->GetOutputPort());
  contour->SetValue(0, dTh1);
  /*
  contour->Update();
  vtkPolyData* polydata = contour->GetOutput();
  polydata->Update();
  bool ret = true;
  if ( polydata->GetNumberOfPoints() < 1 ||
      polydata->GetNumberOfCells() < 1 )
  {
    vtkPolyDataMapper* mapper = vtkPolyDataMapper::SafeDownCast( actor_out->GetMapper() );
    mapper->SetInput( polydata );
    ret = false;
  }
  else*/
  {
    vtkSmartPointer<vtkPolyDataConnectivityFilter> conn = vtkSmartPointer<vtkPolyDataConnectivityFilter>::New();
    conn->SetInputConnection( contour->GetOutputPort() );
    conn->SetExtractionModeToLargestRegion();
    vtkSmartPointer<vtkSmoothPolyDataFilter> smoother = vtkSmartPointer<vtkSmoothPolyDataFilter>::New();
    if ( bAllRegions )
      smoother->SetInputConnection( contour->GetOutputPort() );
    else
      smoother->SetInputConnection( conn->GetOutputPort() );
    smoother->SetNumberOfIterations( nSmoothIterations );
    smoother->FeatureEdgeSmoothingOn();
    smoother->SetEdgeAngle( 90 );
    vtkSmartPointer<vtkPolyDataNormals> normals = vtkSmartPointer<vtkPolyDataNormals>::New();
    normals->SetInputConnection( smoother->GetOutputPort() );
//    normals->SetInput( polydata );
    normals->SetFeatureAngle( 90 );
    vtkSmartPointer<vtkStripper> stripper = vtkSmartPointer<vtkStripper>::New();
    stripper->SetInputConnection( normals->GetOutputPort() );
    vtkPolyDataMapper* mapper = vtkPolyDataMapper::SafeDownCast( actor_out->GetMapper() );
    mapper->SetInputConnection( stripper->GetOutputPort() );
    mapper->ScalarVisibilityOn();
  }

  return true;
}

bool MyVTKUtils::BuildVolume( vtkImageData* data_in,
                           double dTh1, double dTh2,
                           vtkVolume* vol_out )
{
  vtkSmartPointer<vtkPiecewiseFunction> tfun =
    vtkSmartPointer<vtkPiecewiseFunction>::New();
  tfun->AddPoint(dTh1-0.001, 0.0);
  tfun->AddPoint(dTh1, 0.8);
  tfun->AddPoint(dTh2, 1.0);

  vtkSmartPointer<vtkColorTransferFunction> ctfun =
    vtkSmartPointer<vtkColorTransferFunction>::New();
  ctfun->AddRGBPoint( 0.0, 0.0, 0.0, 0.0 );
  ctfun->AddRGBPoint( dTh1, 0.25, 0.25, 0.25 );
  ctfun->AddRGBPoint( (dTh1+dTh2) / 2, 0.4, 0.4, 0.4 );
  ctfun->AddRGBPoint( dTh2, 1, 1, 1 );

  /* vtkSmartPointer<vtkVolumeRayCastCompositeFunction> compositeFunction =
     vtkSmartPointer<vtkVolumeRayCastCompositeFunction>::New();

   vtkSmartPointer<vtkVolumeRayCastMapper> volumeMapper = vtkSmartPointer<vtkVolumeRayCastMapper>::New();
   volumeMapper->SetVolumeRayCastFunction(compositeFunction);
  // vtkOpenGLVolumeShearWarpMapper* volumeMapper = vtkOpenGLVolumeShearWarpMapper::New();
  // vtkOpenGLVolumeTextureMapper3D* volumeMapper = vtkOpenGLVolumeTextureMapper3D::New();
   // vtkVolumeTextureMapper2D* volumeMapper = vtkVolumeTextureMapper2D::New();
  */
  vtkSmartPointer<vtkFixedPointVolumeRayCastMapper> volumeMapper =
    vtkSmartPointer<vtkFixedPointVolumeRayCastMapper>::New();

  vtkSmartPointer<vtkImageCast> cast = vtkSmartPointer<vtkImageCast>::New();
  cast->SetInput( data_in );
  cast->SetOutputScalarTypeToUnsignedShort();

// qDebug() << volumeMapper->GetIntermixIntersectingGeometry();
  volumeMapper->SetInputConnection( cast->GetOutputPort() );
  volumeMapper->SetSampleDistance(0.25);
  volumeMapper->SetMaximumImageSampleDistance(5);
// volumeMapper->SetCroppingRegionPlanes(0, dim[0]*2, 0, dim[1]*2, 0, dim[2]*2-16*2);
// volumeMapper->CroppingOn();

  vtkSmartPointer<vtkVolumeProperty> volumeProperty =
    vtkSmartPointer<vtkVolumeProperty>::New();
  volumeProperty->SetColor(ctfun);
  volumeProperty->SetScalarOpacity(tfun);
  volumeProperty->SetInterpolationTypeToLinear();
  volumeProperty->ShadeOff();

  vol_out->SetMapper( volumeMapper );
  vol_out->SetProperty( volumeProperty );

  return true;
}

void MyVTKUtils::GetLivewirePoints( vtkImageData* image_in,
                                 int nPlane_in, int nSlice_in,
                                 double* pt1_in, double* pt2_in,
                                 vtkPoints* pts_out )
{
  vtkSmartPointer<vtkImageClip> m_imageClip =
    vtkSmartPointer<vtkImageClip>::New();
  vtkSmartPointer<vtkDijkstraImageGeodesicPath> m_path =
    vtkSmartPointer<vtkDijkstraImageGeodesicPath>::New();
  vtkSmartPointer<vtkImageChangeInformation> m_info =
    vtkSmartPointer<vtkImageChangeInformation>::New();
  int m_nPlane = nPlane_in;
  int m_nSlice = nSlice_in;

  m_imageClip->SetInput( image_in );
  int ext[6];
  image_in->GetExtent( ext );
  ext[m_nPlane*2] = ext[m_nPlane*2 + 1] = m_nSlice;
  m_imageClip->SetOutputWholeExtent( ext );
  m_imageClip->ClipDataOn();
  m_imageClip->ReleaseDataFlagOff();
  m_imageClip->Update();

  vtkSmartPointer<vtkImageAnisotropicDiffusion2D> smooth =
    vtkSmartPointer<vtkImageAnisotropicDiffusion2D>::New();
  smooth->SetInputConnection( m_imageClip->GetOutputPort() );
  smooth->SetDiffusionFactor( 0.75 );
  smooth->SetDiffusionThreshold( 50.0 );
  smooth->SetNumberOfIterations( 5 );

  /* vtkSmartPointer<vtkImageGaussianSmooth> smooth =
     vtkSmartPointer<vtkImageGaussianSmooth>::New();
   smooth->SetInputConnection( clip->GetOutputPort() );
   smooth->SetStandardDeviations( 1, 1, 1 );*/

  vtkSmartPointer<vtkImageGradientMagnitude> grad =
    vtkSmartPointer<vtkImageGradientMagnitude>::New();
  grad->SetDimensionality( 2 );
  grad->HandleBoundariesOn();
  grad->SetInputConnection( smooth->GetOutputPort() );
  grad->Update();

  double* range = grad->GetOutput()->GetScalarRange();
  vtkSmartPointer<vtkImageShiftScale> scale =
    vtkSmartPointer<vtkImageShiftScale>::New();
  scale->SetShift( -1.0*range[1] );
  scale->SetScale( 255.0 /( range[0] - range[1] ) );
  scale->SetOutputScalarTypeToShort();
  scale->SetInputConnection( grad->GetOutputPort() );
  scale->ReleaseDataFlagOff();

  m_info->SetInputConnection( scale->GetOutputPort() );
  int n[3] = { 0, 0, 0 };
  n[m_nPlane] = -1*m_nSlice;
  m_info->SetExtentTranslation( n );
  m_info->Update();

  vtkImageData* m_imageSlice = scale->GetOutput();
  m_path->SetInputConnection( m_info->GetOutputPort() );
  // m_path->Update();

  double pt1[3], pt2[3];
// double* orig = image_in->GetOrigin();
  for ( int i = 0; i < 3; i++ )
  {
    // pt1[i] = pt1_in[i] - orig[i];
    // pt2[i] = pt2_in[i] - orig[i];
    pt1[i] = pt1_in[i];
    pt2[i] = pt2_in[i];
  }

  vtkIdType beginVertId = m_imageSlice->FindPoint( pt1 );
  vtkIdType endVertId = m_imageSlice->FindPoint( pt2 );
//  cout << beginVertId << "  " << endVertId << endl;

  if ( beginVertId == -1 || endVertId == -1 )
  {
    // cout << "can not find point: " << pt1_in[0] << " " << pt1_in[1] << " " << pt1_in[2] << ", "
    //   << pt2_in[0] << " " << pt2_in[1] << " " << pt2_in[2] << endl;
    return;
  }

  m_path->SetStartVertex( endVertId );
  m_path->SetEndVertex( beginVertId );
  m_path->Update();

  vtkPolyData *pd = m_path->GetOutput();
  vtkIdType npts = 0, *pts = NULL;
  pd->GetLines()->InitTraversal();
  pd->GetLines()->GetNextCell( npts, pts );
//  cout << npts << endl;
  double offset[3] = { 0, 0, 0 };
  double* vs = image_in->GetSpacing();
  offset[m_nPlane] = m_nSlice*vs[m_nPlane];
  for ( int i = 0; i < npts; i++ )
  {
    double* p = pd->GetPoint( pts[i] );
    // cout << p[0] << " " << p[1] << " " << p[2] << endl;
    pts_out->InsertNextPoint( p[0] + offset[0],
                              p[1] + offset[1],
                              p[2] + offset[2] );
  }

}
