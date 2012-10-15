/**
 * @file  kvlConstructMultiResolutionMeshCollection.cxx
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
#include "kvlMultiResolutionAtlasMesher.h"
#include "kvlCompressionLookupTable.h"
#include "itkImageFileReader.h"
#include <fstream>


int main( int argc, char** argv )
{

  // Sanity check on input
  if ( argc < 8 )
  {
    std::cerr << "Usage: " << argv[ 0 ] << " numberOfUpsamplingSteps meshSizeX meshSizeY meshSizeZ stiffness"
              << " tryToBeSparse fileName1 [ fileName2 ... ]" << std::endl;

    return -1;
  }

  // Retrieve the input parameters
  std::ostringstream  inputParserStream;
  for ( int argumentNumber = 1; argumentNumber < 7; argumentNumber++ )
  {
    inputParserStream << argv[ argumentNumber ] << " ";
  }
  std::istringstream  inputStream( inputParserStream.str().c_str() );
  unsigned int  numberOfUpsamplingSteps;
  unsigned int  meshSizeX;
  unsigned int  meshSizeY;
  unsigned int  meshSizeZ;
  float  stiffness;
  bool  tryToBeSparse;
  inputStream >> numberOfUpsamplingSteps >> meshSizeX >> meshSizeY >> meshSizeZ >> stiffness >> tryToBeSparse;

  try
  {
    // Read the images (in original ushort format)
    typedef kvl::CompressionLookupTable::ImageType  InputImageType;
    std::vector< InputImageType::ConstPointer >  originalImages;
    for ( int argumentNumber = 7; argumentNumber < argc; argumentNumber++ )
    {
      // Read the input image
      typedef itk::ImageFileReader< InputImageType >  ReaderType;
      ReaderType::Pointer  reader = ReaderType::New();
      reader->SetFileName( argv[ argumentNumber ] );
      reader->Update();
      InputImageType::ConstPointer  originalImage = reader->GetOutput();

      // Over-ride the spacing and origin since at this point we can't deal with that
      const double spacing[] = { 1, 1, 1 };
      const double origin[] = { 0, 0, 0 };
      const_cast< InputImageType* >( originalImage.GetPointer() )->SetSpacing( spacing );
      const_cast< InputImageType* >( originalImage.GetPointer() )->SetOrigin( origin );

      // Remember this image
      originalImages.push_back( originalImage );
    }

    // Build a lookup table that maps the original intensities onto uchar starting
    // at 0 and densely packed
    kvl::CompressionLookupTable::Pointer  compressor = kvl::CompressionLookupTable::New();
    compressor->Construct( originalImages );
    compressor->Write( "compressionLookupTable.txt" );

    // Collect the label images resulting from pushing the original images through the
    // lookup table
    typedef kvl::CompressionLookupTable::CompressedImageType  OutputImageType;
    std::vector< OutputImageType::ConstPointer >  labelImages;
    for ( std::vector< InputImageType::ConstPointer >::const_iterator it = originalImages.begin();
          it != originalImages.end(); ++it )
    {
      labelImages.push_back( ( compressor->CompressImage( ( *it ).GetPointer() ) ).GetPointer() );
    }


    // Set up a multi resolution atlas mesher
    kvl::MultiResolutionAtlasMesher::Pointer  mesher = kvl::MultiResolutionAtlasMesher::New();
    mesher->SetLabelImages( labelImages );
    mesher->SetNumberOfUpsamplingSteps( numberOfUpsamplingSteps );
    mesher->SetTryToBeSparse( tryToBeSparse );
    const unsigned int  size[] = { meshSizeX, meshSizeY, meshSizeZ };
    const float  stiffnesses[] = { stiffness, stiffness, stiffness };
    mesher->SetUp( size, stiffnesses );
    mesher->GetCurrentMeshCollection()->Write( "multiResolutionSetup.txt" );
    mesher->Go();
    mesher->GetCurrentMeshCollection()->Write( "multiResolutionDone.txt" );


  }
  catch( itk::ExceptionObject& e )
  {
    std::cerr << e << std::endl;
  }

  return 0;
};

