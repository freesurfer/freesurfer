/**
 * @file  kvlMaskImage.cxx
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
#include "itkMGHImageIOFactory.h"
#include "kvlProgressReporter.h"
#include "itkMaskImageFilter.h"
#include "itkMaskNegatedImageFilter.h"


// Some typedefs
typedef itk::Image< short, 3 >  ImageType;
typedef itk::ImageFileReader< ImageType >  ReaderType;
typedef itk::MaskImageFilter< ImageType, ImageType, ImageType >  MaskerType;
typedef itk::MaskNegatedImageFilter< ImageType, ImageType, ImageType >  NegatedMaskerType;
typedef itk::ImageFileWriter< ImageType >  WriterType;


int main( int argc, char* argv[] )
{
  // Parse input
  if ( argc < 3 )
  {
    std::cerr << argv[ 0 ] << " imageFileName maskFileName [ invertMask=0 ]" << std::endl;
    return -1;
  }
  const std::string  imageFileName( argv[ 1 ] );
  const std::string  maskFileName( argv[ 2 ] );
  bool invertMask = false;
  if ( argc > 3 )
  {
    std::istringstream invertMaskStream( argv[ 3 ] );
    invertMaskStream >> invertMask;
  }

  // Add support for MGH file format to ITK. An alternative way to add this by default would be
  // to edit ITK's itkImageIOFactory.cxx and explicitly adding it in the code there.
  itk::ObjectFactoryBase::RegisterFactory( itk::MGHImageIOFactory::New() );


  try
  {
    // Read image and mask
    ReaderType::Pointer  reader = ReaderType::New();

    reader->SetFileName( imageFileName.c_str() );
    reader->Update();
    ImageType::Pointer  image = reader->GetOutput();
    image->DisconnectPipeline();

    reader->SetFileName( maskFileName.c_str() );
    reader->Update();
    ImageType::ConstPointer  mask = reader->GetOutput();


    // Mask the image
    itk::ImageToImageFilter< ImageType, ImageType >::Pointer  masker = 0;
    if ( !invertMask )
    {
      MaskerType::Pointer  tmpMasker = MaskerType::New();
      tmpMasker->SetInput1( image );
      tmpMasker->SetInput2( mask );
      masker = tmpMasker;
    }
    else
    {
      NegatedMaskerType::Pointer  tmpMasker = NegatedMaskerType::New();
      tmpMasker->SetInput1( image );
      tmpMasker->SetInput2( mask );
      masker = tmpMasker;
    }
    kvl::ProgressReporter  reporter( masker, "Masking" );


    // Write out
    std::ostringstream  outputFileNameStream;
    outputFileNameStream << itksys::SystemTools::GetFilenameWithoutExtension( imageFileName.c_str() )
                         << "_masked.mgz";
    const std::string  outputFileName = outputFileNameStream.str();
    WriterType::Pointer  writer = WriterType::New();
    writer->SetInput( masker->GetOutput() );
    writer->SetFileName( outputFileName.c_str() );
    writer->Update();

    std::cout << "Wrote masked image to file " << outputFileName << std::endl;
  }
  catch ( itk::ExceptionObject& e )
  {
    std::cerr << e << std::endl;
    return -1;
  }

  return 0;
};

