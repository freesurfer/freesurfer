/**
 * @file  kvlAutoCrop.cxx
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
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkMGHImageIOFactory.h"
#include "itkLabelObject.h"
#include "itkLabelMap.h"
#include "itkLabelImageToLabelMapFilter.h"
#include "itkAutoCropLabelMapFilter.h"
#include "itkLabelMapToLabelImageFilter.h"



int main( int argc, char* argv[] )
{
  // Sanity check on input
  if ( argc < 2 )
  {
    std::cerr << argv[0] << " inputImage [ borderWidth=0 ]" << std::endl;
    return -1;
  }

  // Parse input
  const std::string  inputImageFileName = argv[ 1 ];
  int borderWidth = 0;
  if ( argc > 2 )
  {
    std::istringstream  borderWidthStream( argv[ 2 ] );
    borderWidthStream >> borderWidth;
  }

  // Add support for MGH file format to ITK. An alternative way to add this by default would be
  // to edit ITK's itkImageIOFactory.cxx and explicitly adding it in the code there.
  itk::ObjectFactoryBase::RegisterFactory( itk::MGHImageIOFactory::New() );

  try
  {
    // Read image
    typedef itk::Image< short, 3 >  ImageType;
    typedef itk::ImageFileReader< ImageType >  ReaderType;
    ReaderType::Pointer  reader = ReaderType::New();
    reader->SetFileName( inputImageFileName.c_str() );

    // Convert to a label map
    typedef itk::LabelObject< short, 3 >  LabelObjectType;
    typedef itk::LabelMap< LabelObjectType >   LabelMapType;
    typedef itk::LabelImageToLabelMapFilter< ImageType, LabelMapType >  ImageToLabelFilterType;
    ImageToLabelFilterType::Pointer  imageToLabelFilter = ImageToLabelFilterType::New();
    imageToLabelFilter->SetInput( reader->GetOutput() );
    imageToLabelFilter->SetBackgroundValue( 0 );

    // Crop label map
    typedef itk::AutoCropLabelMapFilter< LabelMapType >  CropperType;
    CropperType::Pointer  cropper = CropperType::New();
    cropper->SetInput( imageToLabelFilter->GetOutput() );
    CropperType::SizeType  border = {{ borderWidth, borderWidth, borderWidth }};
    cropper->SetCropBorder( border );


    // Convert to image
    typedef itk::LabelMapToLabelImageFilter< LabelMapType, ImageType >  LabelToImageFilterType;
    LabelToImageFilterType::Pointer  labelToImageFilter = LabelToImageFilterType::New();
    labelToImageFilter->SetInput( cropper->GetOutput() );

    // Write out
    std::ostringstream  outputFileNameStream;
    outputFileNameStream << itksys::SystemTools::GetFilenameWithoutExtension( inputImageFileName.c_str() )
                         << "_autoCropped.mgz";
    const std::string  outputFileName = outputFileNameStream.str();
    typedef itk::ImageFileWriter< ImageType >  WriterType;
    WriterType::Pointer  writer = WriterType::New();
    writer->SetInput( labelToImageFilter->GetOutput() );
    writer->SetFileName( outputFileName.c_str() );
    writer->Update();

    std::cout << "Wrote out " << outputFileName << std::endl;

  }
  catch ( itk::ExceptionObject& e )
  {
    std::cerr << e << std::endl;
    return -1;
  }


  return 0;
};

