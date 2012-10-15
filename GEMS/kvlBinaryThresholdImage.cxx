/**
 * @file  kvlBinaryThresholdImage.cxx
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
#include "kvlRegisterer.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkMGHImageIOFactory.h"
#include "kvlProgressReporter.h"
#include "itkBinaryThresholdImageFilter.h"


typedef itk::Image< unsigned short, 3 >  InputImageType;
typedef itk::Image< unsigned char, 3 >  OutputImageType;
typedef itk::BinaryThresholdImageFilter< InputImageType, OutputImageType >  ThresholderType;
typedef itk::ImageFileReader< InputImageType >  ReaderType;
typedef itk::ImageFileWriter< OutputImageType >  WriterType;




int main( int argc, char* argv[] )
{
  // Sanity check on input
  if ( argc < 3 )
  {
    std::cerr << argv[0] << " imageFileName lowerThreshold [ upperThreshold=lowerThreshold insideValue=255 outsideValue=0 ]" << std::endl;
    return -1;
  }


  // Parse input
  const std::string  imageFileName = argv[ 1 ];
  ThresholderType::InputPixelType lowerThreshold;
  std::istringstream  lowerThresholdStream( argv[ 2 ] );
  lowerThresholdStream >> lowerThreshold;

  ThresholderType::InputPixelType  upperThreshold = lowerThreshold;
  if ( argc > 3 )
  {
    std::istringstream  upperThresholdStream( argv[ 3 ] );
    upperThresholdStream >> upperThreshold;
  }

  ThresholderType::OutputPixelType  insideValue = 255;
  if ( argc > 4 )
  {
    int  insideValueInteger;
    std::istringstream  insideValueStream( argv[ 4 ] );
    insideValueStream >> insideValueInteger;
    insideValue = static_cast< ThresholderType::OutputPixelType >( insideValueInteger );
  }

  ThresholderType::OutputPixelType  outsideValue = 0;
  if ( argc > 5 )
  {
    int  outsideValueInteger;
    std::istringstream  outsideValueStream( argv[ 5 ] );
    outsideValueStream >> outsideValueInteger;
    outsideValue = static_cast< ThresholderType::OutputPixelType >( outsideValueInteger );
  }


  // Print out what we have
  std::cout << "imageFileName: " << imageFileName << std::endl;
  std::cout << "lowerThreshold: " << lowerThreshold << std::endl;
  std::cout << "upperThreshold: " << upperThreshold << std::endl;
  std::cout << "insideValue: " << static_cast< int >( insideValue ) << std::endl;
  std::cout << "outsideValue: " << static_cast< int >( outsideValue ) << std::endl;


  // Add support for MGH file format to ITK. An alternative way to add this by default would be
  // to edit ITK's itkImageIOFactory.cxx and explicitly adding it in the code there.
  itk::ObjectFactoryBase::RegisterFactory( itk::MGHImageIOFactory::New() );


  try
  {
    // Read image
    ReaderType::Pointer  reader = ReaderType::New();
    reader->SetFileName( imageFileName.c_str() );

    // Threshold image
    ThresholderType::Pointer  thresholder = ThresholderType::New();
    thresholder->SetInput( reader->GetOutput() );
    thresholder->SetLowerThreshold( lowerThreshold );
    thresholder->SetUpperThreshold( upperThreshold );
    thresholder->SetOutsideValue( outsideValue );
    thresholder->SetInsideValue( insideValue );

    kvl::ProgressReporter  reporter( thresholder, "Thresholding" );

    // Write image
    std::ostringstream  outputFileNameStream;
    outputFileNameStream << itksys::SystemTools::GetFilenameWithoutExtension( imageFileName.c_str() )
                         << "_thresholded.mgz";
    const std::string  outputFileName = outputFileNameStream.str();
    WriterType::Pointer  writer = WriterType::New();
    writer->SetInput( thresholder->GetOutput() );
    writer->SetFileName( outputFileName.c_str() );

    // Update the pipeline
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

