/**
 * @file  kvlOverwriteEmptySlices.cxx
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: Koen Van Leemput
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2012/10/15 21:17:40 $
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
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"



int main( int argc, char* argv[] )
{

  // Sanity check on input
  if ( argc < 3 )
  {
    std::cerr << "Usage: " << argv[ 0 ] << " imageWithEmptySlices imageToBeOverwritten" << std::endl;
    exit( -1 );
  }

  // Collect the input arguments
  const std::string  imageWithEmptySlicesFileName( argv[ 1 ] );
  const std::string  imageToBeOverwrittenFileName( argv[ 2 ] );

  try
  {
    // Read images
    typedef itk::Image< unsigned short, 3 >  ImageType;
    typedef itk::ImageFileReader< ImageType >  ReaderType;
    ReaderType::Pointer  reader = ReaderType::New();

    reader->SetFileName( imageWithEmptySlicesFileName.c_str() );
    reader->Update();
    ImageType::Pointer  imageWithEmptySlices = reader->GetOutput();
    imageWithEmptySlices->DisconnectPipeline();

    reader->SetFileName( imageToBeOverwrittenFileName.c_str() );
    reader->Update();
    ImageType::Pointer  imageToBeOverwritten = reader->GetOutput();
    imageToBeOverwritten->DisconnectPipeline();


    // Loop over all slices in the imageWithEmptySlices, and determine which onces are empty
    ImageType::IndexType  index = imageWithEmptySlices->GetBufferedRegion().GetIndex();
    ImageType::SizeType  size = imageWithEmptySlices->GetBufferedRegion().GetSize();
    const int  numberOfSlices = size[ 2 ];
    std::vector< bool >  sliceIsEmptyArray( numberOfSlices, true );
    for ( int sliceNumber = 0; sliceNumber < numberOfSlices; sliceNumber++ )
    {
      // Construct region corresponding to this slice
      index[ 2 ] = sliceNumber;
      size[ 2 ] = 1;
      ImageType::RegionType  sliceRegion( index, size );

      for ( itk::ImageRegionConstIterator< ImageType >  it( imageWithEmptySlices, sliceRegion );
            !it.IsAtEnd(); ++it )
      {
        if ( it.Value() != 0 )
        {
          sliceIsEmptyArray[ sliceNumber ] = false;
        }
      } // End loop over all pixels in slice

    } // End loop over all slices


    // Loop over all slices in imageToBeOverwritten, and blank out empty slices
    for ( int sliceNumber = 0; sliceNumber < numberOfSlices; sliceNumber++ )
    {
      if ( sliceIsEmptyArray[ sliceNumber ] == false )
      {
        continue;
      }

      // Construct region corresponding to this slice
      index[ 2 ] = sliceNumber;
      size[ 2 ] = 1;
      ImageType::RegionType  sliceRegion( index, size );

      for ( itk::ImageRegionIterator< ImageType >  it( imageToBeOverwritten, sliceRegion );
            !it.IsAtEnd(); ++it )
      {
        it.Value() = 0;
      } // End loop over all pixels in slice

    } // End loop over all slices


    // Write out result
    typedef itk::ImageFileWriter< ImageType >  WriterType;
    WriterType::Pointer  writer = WriterType::New();
    writer->SetFileName( "overwritten.mhd" );
    writer->SetInput( imageToBeOverwritten );
    writer->Update();

  }
  catch( itk::ExceptionObject& e )
  {
    std::cerr << e << std::endl;
  }

  return 0;
};
