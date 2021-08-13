/**
 * @brief Misc utility class.
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
#include <vtkMarchingCubes.h>
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
#include <vtkVolumeProperty.h>
#include <vtkVolume.h>
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
#include <vtkCleanPolyData.h>
#include <vtkImageResample.h>
#include <vtkImageReslice.h>
#include <vtkWindowedSincPolyDataFilter.h>
#include <QFileInfo>
#include <QDebug>
#include <QMap>
#include "vtkDataSetMapper.h"
#include "vtkGlyph3DMapper.h"
#include "vtkThreshold.h"
#include "vtkDataSetAttributes.h"
#include "vtkGeometryFilter.h"
#include "vtkPointData.h"
#include "vtkFloatArray.h"
#include "vtkPassThrough.h"
#include "vtkDiscreteMarchingCubes.h"
#if VTK_MAJOR_VERSION > 5
#include "vtkFlyingEdges3D.h"
#endif

bool MyVTKUtils::VTKScreenCapture( vtkRenderWindow* renderWnd,
                                   vtkRenderer* renderer,
                                   const char* filename,
                                   bool bAntiAliasing,
                                   int nMag )
{
  Q_UNUSED(bAntiAliasing);
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
  {
    writer = vtkJPEGWriter::New();
  }
  else if ( ext == "bmp" )
  {
    writer = vtkBMPWriter::New();
  }
  else if ( ext == "ps" )
  {
    writer = vtkPostScriptWriter::New();
  }
  else if ( ext == "tif" || ext == "tiff" )
  {
    writer = vtkTIFFWriter::New();
  }
  else
  {
    writer = vtkPNGWriter::New();
    if ( ext != "png" )
    {
      fn += ".png";
    }
  }

  bool ret = true;
  if (writer)
  {
    // bool bCurrentAA = GetAntialiasing() > 0;
    // SetAntialiasing(bAntiAliasing, false);
    vtkRenderLargeImage* image = vtkRenderLargeImage::New();
    image->SetInput( renderer );
    image->SetMagnification( nMag );
#if VTK_MAJOR_VERSION > 5
    writer->SetInputData( image->GetOutput() );
#else
    writer->SetInput( image->GetOutput() );
#endif
    writer->SetFileName( fn.toUtf8().data() );
    writer->Write();
    if ( writer->GetErrorCode() != 0 )
    {
      ret = false;
    }
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

bool MyVTKUtils::BuildLabelContourActor( vtkImageData* data_in,
                                         int labelIndex,
                                         vtkActor* actor_out, int nSmoothIterations, int* ext, bool bAllRegions, bool bUpsample, bool bVoxelized, bool bDilate)
{
  Q_UNUSED(ext);
  int i = labelIndex;
  vtkSmartPointer<vtkImageThreshold> threshold = vtkSmartPointer<vtkImageThreshold>::New();
#if VTK_MAJOR_VERSION > 5
  threshold->SetInputData( data_in );
#else
  threshold->SetInput( data_in );
#endif
  threshold->ThresholdBetween( i-0.5, i+0.5 );
  threshold->ReplaceOutOn();
  threshold->SetOutValue( 0 );
  vtkSmartPointer<vtkImageReslice> resampler = vtkSmartPointer<vtkImageReslice>::New();
  if (bUpsample)
  {
//    resampler->SetAxisMagnificationFactor(0, 2.0);
//    resampler->SetAxisMagnificationFactor(1, 2.0);
//    resampler->SetAxisMagnificationFactor(2, 2.0);
    double vs[3];
    data_in->GetSpacing(vs);
    resampler->SetOutputSpacing(vs[0]/2, vs[1]/2, vs[2]/2);
    resampler->SetInputConnection(threshold->GetOutputPort());
  }

  vtkSmartPointer<vtkMarchingCubes> contour = vtkSmartPointer<vtkMarchingCubes>::New();
  if (bDilate && !bVoxelized)
  {
    int nSwell = 2;
    vtkSmartPointer<vtkImageDilateErode3D> dilate = vtkSmartPointer<vtkImageDilateErode3D>::New();
    dilate->SetInputConnection(bUpsample? resampler->GetOutputPort() : threshold->GetOutputPort());
    dilate->SetKernelSize(nSwell, nSwell, nSwell);
    dilate->SetDilateValue(labelIndex);
    dilate->SetErodeValue(0);
    contour->SetInputConnection(dilate->GetOutputPort());
  }
  else
    contour->SetInputConnection(bUpsample? resampler->GetOutputPort() : threshold->GetOutputPort());

  contour->SetValue(0, i);

  vtkSmartPointer<vtkPolyDataConnectivityFilter> conn = vtkSmartPointer<vtkPolyDataConnectivityFilter>::New();
  conn->SetInputConnection( contour->GetOutputPort() );
  conn->SetExtractionModeToLargestRegion();
  vtkSmartPointer<vtkWindowedSincPolyDataFilter> smoother = vtkSmartPointer<vtkWindowedSincPolyDataFilter>::New();
  if ( bAllRegions )
  {
    smoother->SetInputConnection( contour->GetOutputPort() );
  }
  else
  {
    smoother->SetInputConnection( conn->GetOutputPort() );
  }
  smoother->SetNumberOfIterations( nSmoothIterations );
  //   smoother->FeatureEdgeSmoothingOn();
  //   smoother->SetEdgeAngle( 90 );
  vtkSmartPointer<vtkPolyDataNormals> normals = vtkSmartPointer<vtkPolyDataNormals>::New();
  normals->SetInputConnection( smoother->GetOutputPort() );
  //    normals->SetInput( polydata );
  normals->SetFeatureAngle( 90 );
  vtkSmartPointer<vtkTriangleFilter> stripper = vtkSmartPointer<vtkTriangleFilter>::New();
  stripper->SetInputConnection( normals->GetOutputPort() );
  vtkSmartPointer<vtkCleanPolyData> cleaner = vtkSmartPointer<vtkCleanPolyData>::New();
  cleaner->SetInputConnection(stripper->GetOutputPort());
//  vtkPolyDataMapper* mapper = vtkPolyDataMapper::SafeDownCast( actor_out->GetMapper() );
  if (!bVoxelized)
  {
    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    cleaner->Update();
    mapper->SetInputConnection( cleaner->GetOutputPort() );
    actor_out->SetMapper(mapper);
    mapper->ScalarVisibilityOn();
  }
  else
  {
    threshold->SetOutputScalarTypeToFloat();
    threshold->Update();
    vtkImageData* outputImage = threshold->GetOutput();
    int* dim = outputImage->GetDimensions();
    double* voxel_size = outputImage->GetSpacing();
    double* origin = outputImage->GetOrigin();
    float* p = static_cast<float*>(outputImage->GetScalarPointer());
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkFloatArray> scalars = vtkSmartPointer<vtkFloatArray>::New();
    for (int i = 0; i < dim[0]; i++)
    {
      for (int j = 0; j < dim[1]; j++)
      {
        for (int k = 0; k < dim[2]; k++)
        {
          float val = p[i+j*dim[0]+k*dim[0]*dim[1]];
          if (val > 0)
          {
            points->InsertNextPoint(origin[0]+voxel_size[0]*i, origin[1]+voxel_size[1]*j, origin[2]+voxel_size[2]*k);
            scalars->InsertNextValue(val);
          }
        }
      }
    }
    vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
    polydata->SetPoints(points);
    polydata->GetPointData()->SetScalars(scalars);
    vtkSmartPointer<vtkCubeSource> cube = vtkSmartPointer<vtkCubeSource>::New();
    cube->SetXLength(voxel_size[0]);
    cube->SetYLength(voxel_size[1]);
    cube->SetZLength(voxel_size[2]);
    vtkSmartPointer<vtkGlyph3DMapper> glyph = vtkSmartPointer<vtkGlyph3DMapper>::New();
    glyph->SetSourceConnection(cube->GetOutputPort());
    glyph->SetScaleModeToNoDataScaling();
#if VTK_MAJOR_VERSION > 5
    glyph->SetInputData(polydata);
#else
    vtkSmartPointer<vtkPassThrough> pass = vtkSmartPointer<vtkPassThrough>::New();
    pass->SetInput(polydata);
    glyph->SetInputConnection(pass->GetOutputPort());
#endif
    glyph->Update();
    actor_out->SetMapper(glyph);
    glyph->ScalarVisibilityOn();
  }

  return true;
}

// test multiple contours
bool MyVTKUtils::BuildLabelContourActor( vtkImageData* data_in,
                                         const QList<int>& labelIndices,
                                         vtkActor* actor_out, int nSmoothIterations, int* ext, bool bAllRegions, bool bUpsample )
{
  Q_UNUSED(ext);
//  double nValue = 1;
//  int nSwell = 2;

  vtkSmartPointer<vtkAppendPolyData> append = vtkSmartPointer<vtkAppendPolyData>::New();
  foreach (int i, labelIndices)
  {
    vtkSmartPointer<vtkImageThreshold> threshold = vtkSmartPointer<vtkImageThreshold>::New();
#if VTK_MAJOR_VERSION > 5
    threshold->SetInputData( data_in );
#else
    threshold->SetInput( data_in );
#endif
    threshold->ThresholdBetween( i-0.5, i+0.5 );
    threshold->ReplaceOutOn();
    threshold->SetOutValue( 0 );
    vtkSmartPointer<vtkImageReslice> resampler = vtkSmartPointer<vtkImageReslice>::New();
    if (bUpsample)
    {
//      resampler->SetAxisMagnificationFactor(0, 2.0);
//      resampler->SetAxisMagnificationFactor(1, 2.0);
//      resampler->SetAxisMagnificationFactor(2, 2.0);
      double vs[3];
      data_in->GetSpacing(vs);
      resampler->SetOutputSpacing(vs[0]/2, vs[1]/2, vs[2]/2);
      resampler->SetInputConnection(threshold->GetOutputPort());
    }
    vtkSmartPointer<vtkMarchingCubes> contour = vtkSmartPointer<vtkMarchingCubes>::New();
    contour->SetInputConnection( bUpsample? resampler->GetOutputPort() : threshold->GetOutputPort());
    contour->SetValue(0, i);
    append->AddInputConnection(contour->GetOutputPort());
  }
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
    conn->SetInputConnection( append->GetOutputPort() );
    conn->SetExtractionModeToLargestRegion();
    vtkSmartPointer<vtkWindowedSincPolyDataFilter> smoother = vtkSmartPointer<vtkWindowedSincPolyDataFilter>::New();
    if ( bAllRegions )
    {
      smoother->SetInputConnection( append->GetOutputPort() );
    }
    else
    {
      smoother->SetInputConnection( conn->GetOutputPort() );
    }
    smoother->SetNumberOfIterations( nSmoothIterations );
    //   smoother->FeatureEdgeSmoothingOn();
    //   smoother->SetEdgeAngle( 90 );
    vtkSmartPointer<vtkPolyDataNormals> normals = vtkSmartPointer<vtkPolyDataNormals>::New();
    normals->SetInputConnection( smoother->GetOutputPort() );
    //    normals->SetInput( polydata );
    normals->SetFeatureAngle( 90 );
    vtkSmartPointer<vtkTriangleFilter> stripper = vtkSmartPointer<vtkTriangleFilter>::New();
    stripper->SetInputConnection( normals->GetOutputPort() );
    vtkSmartPointer<vtkCleanPolyData> cleaner = vtkSmartPointer<vtkCleanPolyData>::New();
    cleaner->SetInputConnection(stripper->GetOutputPort());
    cleaner->Update();
    vtkPolyDataMapper* mapper = vtkPolyDataMapper::SafeDownCast( actor_out->GetMapper() );
    mapper->SetInputConnection( cleaner->GetOutputPort() );
    mapper->ScalarVisibilityOn();
  }

  return true;
}

bool MyVTKUtils::BuildContourActor( vtkImageData* data_in,
                                    double dTh1, double dTh2,
                                    vtkActor* actor_out, int nSmoothIterations, int* ext, bool bAllRegions,
                                    bool bUpsample, bool bDilate)
{
  double nValue = 1;
  vtkSmartPointer<vtkImageThreshold> threshold = vtkSmartPointer<vtkImageThreshold>::New();

  vtkSmartPointer<vtkImageReslice> resampler = vtkSmartPointer<vtkImageReslice>::New();
//  resampler->SetAxisMagnificationFactor(0, 2.0);
//  resampler->SetAxisMagnificationFactor(1, 2.0);
//  resampler->SetAxisMagnificationFactor(2, 2.0);
  double vs[3];
  data_in->GetSpacing(vs);
  resampler->SetOutputSpacing(vs[0]/2, vs[1]/2, vs[2]/2);
  if ( ext )
  {
    vtkSmartPointer<vtkImageClip> clipper = vtkSmartPointer<vtkImageClip>::New();
    if (bUpsample)
    {
#if VTK_MAJOR_VERSION > 5
      resampler->SetInputData(data_in);
#else
      resampler->SetInput(data_in);
#endif
      clipper->SetInputConnection(resampler->GetOutputPort());
    }
    else
#if VTK_MAJOR_VERSION > 5
      clipper->SetInputData( data_in );
#else
      clipper->SetInput( data_in );
#endif
    clipper->SetOutputWholeExtent( ext );
    threshold->SetInputConnection( clipper->GetOutputPort() );
  }
  else
  {
    if (bUpsample)
    {
#if VTK_MAJOR_VERSION > 5
      resampler->SetInputData(data_in);
#else
      resampler->SetInput(data_in);
#endif
      threshold->SetInputConnection(resampler->GetOutputPort());
    }
    else
#if VTK_MAJOR_VERSION > 5
      threshold->SetInputData( data_in );
#else
      threshold->SetInput( data_in );
#endif
  }
  threshold->ThresholdByLower( dTh2 );
  threshold->ReplaceOutOn();
  threshold->SetOutValue( dTh1-0.00001 );

  /*
  vtkSmartPointer<vtkImageDilateErode3D> erode = vtkSmartPointer<vtkImageDilateErode3D>::New();
  erode->SetInputConnection(dilate->GetOutputPort());
  erode->SetKernelSize(1, 1, 1);
  erode->SetDilateValue(0);
  erode->SetErodeValue(nValue);
  */

  vtkSmartPointer<vtkContourFilter> contour = vtkSmartPointer<vtkContourFilter>::New();
  if (bDilate)
  {
    int nSwell = 2;
    vtkSmartPointer<vtkImageDilateErode3D> dilate = vtkSmartPointer<vtkImageDilateErode3D>::New();
    dilate->SetInputConnection(threshold->GetOutputPort());
    dilate->SetKernelSize(nSwell, nSwell, nSwell);
    dilate->SetDilateValue(dTh2);
    dilate->SetErodeValue(0);
    contour->SetInputConnection( dilate->GetOutputPort());
  }
  else
    contour->SetInputConnection( threshold->GetOutputPort());
  contour->SetValue(0, dTh1);
  {
    vtkSmartPointer<vtkPolyDataConnectivityFilter> conn = vtkSmartPointer<vtkPolyDataConnectivityFilter>::New();
    conn->SetInputConnection( contour->GetOutputPort() );
    conn->SetExtractionModeToLargestRegion();
    vtkSmartPointer<vtkWindowedSincPolyDataFilter> smoother = vtkSmartPointer<vtkWindowedSincPolyDataFilter>::New();
    if ( bAllRegions )
    {
      smoother->SetInputConnection( contour->GetOutputPort() );
    }
    else
    {
      smoother->SetInputConnection( conn->GetOutputPort() );
    }
    smoother->SetNumberOfIterations( nSmoothIterations );
    //   smoother->SetRelaxationFactor(smoother->GetRelaxationFactor()*2);
    //   smoother->FeatureEdgeSmoothingOn();
    //   smoother->SetEdgeAngle( 90 );
    vtkSmartPointer<vtkPolyDataNormals> normals = vtkSmartPointer<vtkPolyDataNormals>::New();
    normals->SetInputConnection( smoother->GetOutputPort() );
    //    normals->SetInput( polydata );
    normals->SetFeatureAngle( 90 );
    vtkSmartPointer<vtkTriangleFilter> stripper = vtkSmartPointer<vtkTriangleFilter>::New();
    stripper->SetInputConnection( normals->GetOutputPort() );
    vtkSmartPointer<vtkCleanPolyData> cleaner = vtkSmartPointer<vtkCleanPolyData>::New();
    cleaner->SetInputConnection(stripper->GetOutputPort());
    cleaner->Update();
    vtkPolyDataMapper* mapper = vtkPolyDataMapper::SafeDownCast( actor_out->GetMapper() );
    mapper->SetInputConnection( cleaner->GetOutputPort() );
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
#if VTK_MAJOR_VERSION > 5
  cast->SetInputData( data_in );
#else
  cast->SetInput( data_in );
#endif
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

#if VTK_MAJOR_VERSION > 5
  m_imageClip->SetInputData( image_in );
#else
  m_imageClip->SetInput( image_in );
#endif
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


double MyVTKUtils::GetImageDataComponent(char* ptr, int* dim, size_t nNumberOfFrames, size_t i, size_t j, size_t k, size_t nframe, int data_type)
{
  switch (data_type)
  {
  case VTK_UNSIGNED_CHAR:
    return ((unsigned char*)ptr)[(k*dim[0]*dim[1]+j*dim[0]+i)*nNumberOfFrames + nframe];
  case VTK_INT:
    return ((int*)ptr)[(k*dim[0]*dim[1]+j*dim[0]+i)*nNumberOfFrames + nframe];
  case VTK_LONG:
    return ((long*)ptr)[(k*dim[0]*dim[1]+j*dim[0]+i)*nNumberOfFrames + nframe];
  case VTK_FLOAT:
    return ((float*)ptr)[(k*dim[0]*dim[1]+j*dim[0]+i)*nNumberOfFrames + nframe];
  case VTK_SHORT:
    return ((short*)ptr)[(k*dim[0]*dim[1]+j*dim[0]+i)*nNumberOfFrames + nframe];
  }
  return 0;
}

void MyVTKUtils::SetImageDataComponent(char* ptr, int* dim, size_t nNumberOfFrames, size_t i, size_t j, size_t k, size_t nframe, int data_type, double val)
{
  switch (data_type)
  {
  case VTK_UNSIGNED_CHAR:
    ((unsigned char*)ptr)[(k*dim[0]*dim[1]+j*dim[0]+i)*nNumberOfFrames + nframe] = ((unsigned char)val);
    break;
  case VTK_INT:
    ((int*)ptr)[(k*dim[0]*dim[1]+j*dim[0]+i)*nNumberOfFrames + nframe] = ((int)val);
    break;
  case VTK_LONG:
    ((long*)ptr)[(k*dim[0]*dim[1]+j*dim[0]+i)*nNumberOfFrames + nframe] = (long)val;
    break;
  case VTK_FLOAT:
    ((float*)ptr)[(k*dim[0]*dim[1]+j*dim[0]+i)*nNumberOfFrames + nframe] = (float)val;
    break;
  case VTK_SHORT:
    ((short*)ptr)[(k*dim[0]*dim[1]+j*dim[0]+i)*nNumberOfFrames + nframe] = (short)val;
    break;
  }
}
