/**
 * @file  kvlCompressImage.cxx
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
#include "kvlCompressionLookupTable.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkMGHImageIOFactory.h"



int main( int argc, char** argv )
{
  // Sanity check on input
  if ( argc < 2 )
  {
    std::cerr << "Usage: "<< argv[ 0 ] << " imageFileName [ lookupFileName ]" << std::endl;
    return -1;
  }

  // Add support for MGH file format to ITK. An alternative way to add this by default would be
  // to edit ITK's itkImageIOFactory.cxx and explicitly adding it in the code there.
  itk::ObjectFactoryBase::RegisterFactory( itk::MGHImageIOFactory::New() );

  // Read the image
  typedef kvl::CompressionLookupTable::ImageType  ImageType;
  typedef itk::ImageFileReader< ImageType >  ReaderType;
  ReaderType::Pointer  reader = ReaderType::New();
  reader->SetFileName( argv[ 1 ] );
  reader->Update();

  // Construct the compressor, either by reading it from file or from the image itself
  kvl::CompressionLookupTable::Pointer  compressor = kvl::CompressionLookupTable::New();
  if ( argc > 2 )
  {
    compressor->Read( argv[ 2 ] );
  }
  else
  {
    std::vector< ImageType::ConstPointer >  images;
    images.push_back( reader->GetOutput() );
    compressor->Construct( images );
    compressor->Write( "compressionLookupTable.txt" );
  }

  // Push the image through the compressor
  typedef kvl::CompressionLookupTable::CompressedImageType  CompressedImageType;
  CompressedImageType::Pointer  compressedImage =  compressor->CompressImage( reader->GetOutput() );

  // Write the compressed image out
  typedef itk::ImageFileWriter< CompressedImageType >  WriterType;
  WriterType::Pointer  writer = WriterType::New();
  writer->SetFileName( "compressed.mhd" );
  writer->SetInput( compressedImage );
  writer->Write();

  return 0;
};
