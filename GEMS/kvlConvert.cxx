/**
 * @file  kvlConvert.cxx
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
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkMGHImageIOFactory.h"
#include "itkImageIOFactory.h"


int main( int argc, char* argv[] )
{
  // Sanity check on input
  if ( argc < 3 )
  {
    std::cerr << argv[0] << " inputImageFileName outputImageFileName [ headerImageFileName ]" << std::endl;
    return -1;
  }

  // Parse input
  const std::string  inputImageFileName( argv[ 1 ] );
  const std::string  outputImageFileName( argv[ 2 ] );
  std::string  headerImageFileName;
  if ( argc > 3 )
  {
    headerImageFileName = argv[ 3 ];
  }


  std::cout << "inputImageFileName: " << inputImageFileName << std::endl;
  std::cout << "outputImageFileName: " << outputImageFileName << std::endl;
  std::cout << "headerImageFileName: " << headerImageFileName << std::endl;


  // Add support for MGH file format to ITK. An alternative way to add this by default would be
  // to edit ITK's itkImageIOFactory.cxx and explicitly adding it in the code there.
  itk::ObjectFactoryBase::RegisterFactory( itk::MGHImageIOFactory::New() );


  try
  {
    // Read the image
    typedef itk::Image< short, 3 >  ImageType;
    typedef itk::ImageFileReader< ImageType >  ReaderType;
    ReaderType::Pointer  reader = ReaderType::New();
    reader->SetFileName( inputImageFileName.c_str() );
    reader->Update();
    ImageType::Pointer  image = reader->GetOutput();

    // Read image to copy header from
    if ( headerImageFileName.size() != 0 )
    {
      // Reader header
      itk::ImageIOBase::Pointer  imageIO = itk::ImageIOFactory::CreateImageIO( headerImageFileName.c_str(), itk::ImageIOFactory::ReadMode );

      if ( imageIO.IsNull() )
      {
        std::cerr << "Could read header of file " << headerImageFileName << std::endl;
        exit( - 1 );
      }
      imageIO->SetFileName( headerImageFileName.c_str() );
      imageIO->ReadImageInformation();


      // Overwrite header info
      double  spacing[ 3 ];
      double origin[ 3 ];
      ImageType::DirectionType  direction;
      for ( int i = 0; i < 3; i++ )
      {
        spacing[ i ] = imageIO->GetSpacing( i );
        origin[ i ] = imageIO->GetOrigin( i );

        const std::vector< double >  axis = imageIO->GetDirection( i );
        for ( int j = 0; j< 3; j++ )
        {
          direction[ j ][ i ] = axis[ j ];
        }
      }
      image->SetSpacing( spacing );
      image->SetOrigin( origin );
      image->SetDirection( direction );

    } // End test if we should read another image to copy header info from


    // Write out result
    typedef itk::ImageFileWriter< ImageType >  WriterType;
    WriterType::Pointer  writer = WriterType::New();
    writer->SetFileName( outputImageFileName.c_str() );
    writer->SetInput( image );
    writer->Update();
    std::cout << "Wrote " << outputImageFileName << std::endl;

  }
  catch ( itk::ExceptionObject& e )
  {
    std::cerr << e << std::endl;
    return -1;
  }


  return 0;
};


