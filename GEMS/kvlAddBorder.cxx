/**
 * @file  kvlAddBorder.cxx
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: Koen Van Leemput
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2012/10/15 21:17:38 $
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
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"



int main( int argc, char* argv[] )
{
  if ( argc < 5 )
  {
    std::cerr << argv[0] << " xBorderWidth yBorderWidth zBorderWith imageFileName" << std::endl;
    return -1;
  }

  // Retrieve the input parameters
  std::ostringstream  inputParserStream;
  for ( int argumentNumber = 1; argumentNumber < 5; argumentNumber++ )
  {
    inputParserStream << argv[ argumentNumber ] << " ";
  }
  std::istringstream  inputStream( inputParserStream.str().c_str() );
  unsigned int  xBorderWidth;
  unsigned int  yBorderWidth;
  unsigned int  zBorderWidth;
  std::string  imageFileName;
  inputStream >> xBorderWidth >> yBorderWidth >> zBorderWidth >> imageFileName;


  // Read the image
  typedef itk::Image< short, 3 >  ImageType;
  typedef itk::ImageFileReader< ImageType >  ReaderType;
  ReaderType::Pointer  reader = ReaderType::New();
  reader->SetFileName( imageFileName.c_str() );
  reader->Update();
  ImageType::ConstPointer  inputImage = reader->GetOutput();


  // Gets the image regions
  ImageType::RegionType  inputImageRegion = inputImage->GetLargestPossibleRegion();

  ImageType::IndexType  outputImageIndex = {{ 0, 0, 0 }};
  ImageType::SizeType  outputImageSize;
  outputImageSize[ 0 ] = inputImageRegion.GetSize()[ 0 ] + 2 * xBorderWidth;
  outputImageSize[ 1 ] = inputImageRegion.GetSize()[ 1 ] + 2 * yBorderWidth;
  outputImageSize[ 2 ] = inputImageRegion.GetSize()[ 2 ] + 2 * zBorderWidth;
  ImageType::RegionType  outputImageRegion( outputImageIndex, outputImageSize );

  ImageType::IndexType  walkingOutputImageIndex = {{ xBorderWidth, yBorderWidth, zBorderWidth }};
  ImageType::SizeType  walkingOutputImageSize;
  walkingOutputImageSize[ 0 ] = inputImageRegion.GetSize()[ 0 ];
  walkingOutputImageSize[ 1 ] = inputImageRegion.GetSize()[ 1 ];
  walkingOutputImageSize[ 2 ] = inputImageRegion.GetSize()[ 2 ];
  ImageType::RegionType  walkingOutputImageRegion( walkingOutputImageIndex, walkingOutputImageSize );


  // Create an empty image, and fill it up
  const double spacing[] = { 1, 1, 1 };
  const double origin[] = { 0, 0, 0 };
  ImageType::Pointer  outputImage = ImageType::New();
  outputImage->SetRegions( outputImageRegion );
  outputImage->SetSpacing( spacing );
  outputImage->SetOrigin( origin );
  outputImage->Allocate();
  outputImage->FillBuffer( 0 );

  itk::ImageRegionConstIterator< ImageType >  inputImageIterator( inputImage, inputImageRegion );
  itk::ImageRegionIterator< ImageType >  outputImageIterator( outputImage, walkingOutputImageRegion );
  for ( ; !inputImageIterator.IsAtEnd(); ++inputImageIterator, ++outputImageIterator )
  {
    outputImageIterator.Value() = inputImageIterator.Value();
  }


  // Write out result
  const std::string  outputFileName = "bordered.mhd";
  typedef itk::ImageFileWriter< ImageType >  WriterType;
  WriterType::Pointer  writer = WriterType::New();
  writer->SetInput( outputImage );
  writer->SetFileName( outputFileName.c_str() );
  writer->Write();

  std::cout << "Wrote result to " << outputFileName << std::endl;

  return 0;

};


