/**
 * @file  kvlCrop.cxx
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
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"


int main( int argc, char** argv )
{
  // Sanity check on input
  if ( argc < 3 )
  {
    std::cerr << "Usage: "<< argv[ 0 ] << " marginFraction inputImage1 [inputImage2 ...]" << std::endl;
    return -1;
  }


  // Read margin
  std::istringstream  marginFractionStream( argv[ 1 ] );
  float  marginFraction;
  marginFractionStream >> marginFraction;
  std::cout << "marginFraction: " << marginFraction << std::endl;

  // Read images
  typedef itk::Image< unsigned short, 3 >  ImageType;
  std::vector< ImageType::ConstPointer >  inputImages;
  for ( int inputNumber = 2; inputNumber < argc; inputNumber++ )
  {
    typedef itk::ImageFileReader< ImageType >  ReaderType;
    ReaderType::Pointer  reader = ReaderType::New();
    reader->SetFileName( argv[ inputNumber ] );
    reader->Update();
    inputImages.push_back( reader->GetOutput() );
  }



  // Loop over all images to find the bounding box with content
  int  min[ 3 ];
  for ( int i = 0; i < 3; i++ )
  {
    min[ i ] = itk::NumericTraits< int >::max();
  }
  int  max[ 3 ];
  for ( int i = 0; i < 3; i++ )
  {
    max[ i ] = itk::NumericTraits< int >::min();
  }

  for ( std::vector< ImageType::ConstPointer >::const_iterator  imageIt = inputImages.begin();
        imageIt != inputImages.end(); imageIt++ )
  {
    for ( itk::ImageRegionConstIteratorWithIndex< ImageType >  it( *imageIt, ( *imageIt )->GetBufferedRegion() );
          !it.IsAtEnd(); ++it )
    {
      if ( it.Value() == 0 )
      {
        continue;
      }

      for ( int i = 0; i < 3; i++ )
      {
        if ( it.GetIndex()[ i ] < min[ i ] )
        {
          min[ i ] = it.GetIndex()[ i ];
        }

        if ( it.GetIndex()[ i ] > max[ i ] )
        {
          max[ i ] = it.GetIndex()[ i ];
        }
      }

    }
  }

  std::cout << "Found min [" << min[ 0 ] << "  " << min[ 1 ] << "  " << min[ 2 ] << "]" << std::endl;
  std::cout << "Found max [" << max[ 0 ] << "  " << max[ 1 ] << "  " << max[ 2 ] << "]" << std::endl;


  // Add some margin as a fraction of the found bounding box
  for ( int i = 0; i < 3; i++ )
  {
    const int  margin = static_cast< int >( marginFraction * ( max[ i ] - min[ i ] ) );

    min[ i ] -= margin;
    max[ i ] += margin;
  }


  // Calculate clipped bounding box
  int  clippedMin[ 3 ];
  int  clippedMax[ 3 ];
  for ( int i = 0; i < 3; i++ )
  {
    if ( min[ i ] < 0 )
    {
      clippedMin[ i ] = 0;
    }
    else
    {
      clippedMin[ i ] = min[ i ];
    }

    if ( max[ i ] > static_cast< int >( inputImages[ 0 ]->GetBufferedRegion().GetSize()[i] - 1 ) )
    {
      clippedMax[ i ] = ( inputImages[ 0 ]->GetBufferedRegion().GetSize()[i] - 1 );
    }
    else
    {
      clippedMax[ i ] = max[ i ];
    }
  }


  std::cout << "Cropping at min [" << min[ 0 ] << "  " << min[ 1 ] << "  " << min[ 2 ] << "]" << std::endl;
  std::cout << "     clippedMin [" << clippedMin[ 0 ] << "  " << clippedMin[ 1 ] << "  " << clippedMin[ 2 ] << "]" << std::endl;
  std::cout << "Cropping at max [" << max[ 0 ] << "  " << max[ 1 ] << "  " << max[ 2 ] << "]" << std::endl;
  std::cout << "     clippedMax [" << clippedMax[ 0 ] << "  " << clippedMax[ 1 ] << "  " << clippedMax[ 2 ] << "]" << std::endl;

  // Gets the image regions
  ImageType::IndexType  inputImageIndex = {{ clippedMin[ 0 ], clippedMin[ 1 ], clippedMin[ 2 ] }};
  ImageType::SizeType  inputImageSize;
  for ( int i = 0; i < 3; i++ )
  {
    inputImageSize[ i ] = clippedMax[ i ] - clippedMin[ i ] + 1;
  }
  ImageType::RegionType  inputImageRegion( inputImageIndex, inputImageSize );

  ImageType::IndexType  nonClippedOutputImageIndex = {{ 0, 0, 0 }};
  ImageType::SizeType  nonClippedOutputImageSize;
  for ( int i = 0; i < 3; i++ )
  {
    nonClippedOutputImageSize[ i ] = max[ i ] - min[ i ] + 1;
  }
  ImageType::RegionType  nonClippedOutputImageRegion( nonClippedOutputImageIndex, nonClippedOutputImageSize );

  ImageType::IndexType  outputImageIndex;
  for ( int i = 0; i < 3; i++ )
  {
    outputImageIndex[ i ] = clippedMin[ i ] - min[ i ];
  }
  ImageType::SizeType  outputImageSize = inputImageSize;
  ImageType::RegionType  outputImageRegion( outputImageIndex, outputImageSize );


  // Loop over all images, and crop each time
  for ( unsigned int imageNumber = 0; imageNumber < inputImages.size(); imageNumber++ )
  {
    // Create an empty image, and fill it up
    ImageType::Pointer  outputImage = ImageType::New();
    outputImage->SetRegions( nonClippedOutputImageRegion );
    outputImage->Allocate();
    outputImage->FillBuffer( 0 );

    itk::ImageRegionConstIterator< ImageType >  inputImageIterator( inputImages[ imageNumber ], inputImageRegion );
    itk::ImageRegionIterator< ImageType >  outputImageIterator( outputImage, outputImageRegion );
    for ( ; !inputImageIterator.IsAtEnd(); ++inputImageIterator, ++outputImageIterator )
    {
      outputImageIterator.Value() = inputImageIterator.Value();
    }


    // Construct output filename
    std::ostringstream  outputfileNameStream;
    outputfileNameStream << "cropped";
    if ( imageNumber < 10 )
    {
      outputfileNameStream << "0";
    }
    outputfileNameStream << imageNumber << ".mhd";
    const std::string  outputFileName = outputfileNameStream.str();

    // Write out result
    typedef itk::ImageFileWriter< ImageType >  WriterType;
    WriterType::Pointer  writer = WriterType::New();
    writer->SetInput( outputImage );
    writer->SetFileName( outputFileName.c_str() );
    writer->Write();
    std::cout << "Wrote to " << outputFileName << std::endl;
  }

  return 0;
};
