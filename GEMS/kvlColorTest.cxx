/**
 * @file  kvlColorTest.cxx
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: Koen Van Leemput
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2012/10/15 21:17:39 $
 *    $Revision: 1.3 $
 *
 * Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
 *
 * Terms and conditions for use, reproduction, distribution and contribution
 * are found in the 'FreeSurfer Software License Agreement' contained
 * in the file 'LICENSE' found in the FreeSurfer distribution, and here:
 *
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
 *
 * Reporting: freesurfer@nmr.mgh.harvard.edu
 *
 */
#include "itkVector.h"
#include "itkRGBAPixel.h"
#include "itkImage.h"
#include "itkVTKImageExport.h"
#include "vtkImageImport.h"
#include "vtkImageData.h"
#include "vtkImageWriter.h"
#include "itkImageFileWriter.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkImageActor.h"
#include "vtkCamera.h"




template <typename ITK_Exporter, typename VTK_Importer>
void ConnectPipelines(ITK_Exporter exporter, VTK_Importer* importer)
{
  importer->SetUpdateInformationCallback(exporter->GetUpdateInformationCallback());
  importer->SetPipelineModifiedCallback(exporter->GetPipelineModifiedCallback());
  importer->SetWholeExtentCallback(exporter->GetWholeExtentCallback());
  importer->SetSpacingCallback(exporter->GetSpacingCallback());
  importer->SetOriginCallback(exporter->GetOriginCallback());
  importer->SetScalarTypeCallback(exporter->GetScalarTypeCallback());
  importer->SetNumberOfComponentsCallback(exporter->GetNumberOfComponentsCallback());
  importer->SetPropagateUpdateExtentCallback(exporter->GetPropagateUpdateExtentCallback());
  importer->SetUpdateDataCallback(exporter->GetUpdateDataCallback());
  importer->SetDataExtentCallback(exporter->GetDataExtentCallback());
  importer->SetBufferPointerCallback(exporter->GetBufferPointerCallback());
  importer->SetCallbackUserData(exporter->GetCallbackUserData());
}


int main( int argc, char* argv[] )
{

  //
  //typedef itk::Vector< unsigned char, 4 >  PixelType;
  typedef itk::RGBAPixel< unsigned char >  PixelType;
  typedef itk::Image< PixelType, 3 >  ImageType;
  ImageType::SizeType  size = {{ 100, 100, 100 }};
  ImageType::Pointer  image = ImageType::New();
  image->SetRegions( size );
  image->Allocate();
  PixelType  defaultPixelValue;
  defaultPixelValue[ 0 ] = 200;
  defaultPixelValue[ 1 ] = 50;
  defaultPixelValue[ 2 ] = 50;
  defaultPixelValue[ 3 ] = 80;
  image->FillBuffer( defaultPixelValue );

  //
  {
    typedef itk::ImageFileWriter< ImageType >  WriterType;
    WriterType::Pointer  writer = WriterType::New();
    //writer->SetFileTypeToASCII();
    writer->SetInput( image );
    writer->SetFileName( "haha.vtk" );
    writer->Write();
  }

  // Connect to VTK
  typedef itk::VTKImageExport< ImageType > ExportFilterType;
  ExportFilterType::Pointer itkExporter = ExportFilterType::New();
  itkExporter->SetInput( image );
  vtkImageImport*  vtkImporter = vtkImageImport::New();
  ConnectPipelines(itkExporter, vtkImporter);

  // VTK stuff
  vtkImporter->Update();
  vtkImporter->GetOutput()->Print( std::cout );
  vtkImageWriter*  writer = vtkImageWriter::New();
  writer->SetInput( vtkImporter->GetOutput() );
  writer->SetFileName( "test.vtk" );
  writer->Write();


  vtkRenderer* aRenderer = vtkRenderer::New();
  vtkRenderWindow* renWin = vtkRenderWindow::New();
  renWin->AddRenderer( aRenderer );
  vtkRenderWindowInteractor* iren = vtkRenderWindowInteractor::New();
  iren->SetRenderWindow( renWin );

  vtkImageActor* coronal = vtkImageActor::New();
  coronal->SetInput( vtkImporter->GetOutput() );
  coronal->SetDisplayExtent(0, 63, 32, 32, 0, 92 );

  vtkCamera* aCamera = vtkCamera::New();
  aCamera->SetViewUp( 0, 0, -1 );
  aCamera->SetPosition( 0, 1, 0 );
  aCamera->SetFocalPoint( 0, 0, 0 );
  aCamera->ComputeViewPlaneNormal();

  aRenderer->AddActor( coronal );

  aRenderer->SetActiveCamera( aCamera );
  aRenderer->ResetCamera();
  aCamera->Dolly( 1.5 );

  aRenderer->SetBackground( 1, 1, 1 );
  renWin->SetSize(640, 480);

  aRenderer->ResetCameraClippingRange();

  iren->Initialize();
  iren->Start();

  return 0;
};
