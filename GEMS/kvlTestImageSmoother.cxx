/**
 * @file  kvlTestImageSmoother.cxx
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
#include "kvlImageSmoother.h"
#include "itkGaussianImageSource.h"
#include "itkImageFileWriter.h"



int main( int argc, char* argv[] )
{

  try
  {
    // Generate test image
    typedef kvl::ImageSmoother::ImageType  ImageType;
    typedef itk::GaussianImageSource< ImageType >  SourceType;
    SourceType::Pointer  source = SourceType::New();
    SourceType::SizeType  size = {{ 64, 64, 64 }};
    source->SetSize( size );
    source->Update();
    ImageType::Pointer  image = source->GetOutput();

    // Generate a mask image
    typedef kvl::ImageSmoother::MaskImageType  MaskImageType;
    MaskImageType::Pointer  maskImage = MaskImageType::New();
    maskImage->SetRegions( image->GetBufferedRegion() );
    maskImage->Allocate();
    maskImage->FillBuffer( true );

    // Write out test image and mask image
    typedef itk::ImageFileWriter< ImageType >  WriterType;
    WriterType::Pointer  writer = WriterType::New();
    writer->SetInput( image );
    writer->SetFileName( "image.mhd" );
    writer->Write();

    typedef itk::ImageFileWriter< MaskImageType >  MaskWriterType;
    MaskWriterType::Pointer  maskWriter = MaskWriterType::New();
    maskWriter->SetInput( maskImage );
    maskWriter->SetFileName( "maskImage.mhd" );

    // Set up the image smoother
    kvl::ImageSmoother::Pointer  smoother = kvl::ImageSmoother::New();
    smoother->SetMaskImage( maskImage );
    smoother->SetImage( image );

    // Now smoth for different order polynomials
    for ( int polynomialOrder = 0; polynomialOrder <= 4; polynomialOrder++ )
    {
      //
      std::cout << "Fitting for polynomialOrder: " << polynomialOrder << std::endl;

      smoother->SetPolynomialOrder( polynomialOrder );
      ImageType::Pointer  smoothedImage = smoother->GetSmoothedImage();

      // Write out
      std::ostringstream  fileNameStream;
      fileNameStream << "smoothedForPolymialOrder_" << polynomialOrder << ".mhd";
      writer->SetInput( smoothedImage );
      writer->SetFileName( fileNameStream.str().c_str() );
      writer->Write();
    }

  }
  catch ( itk::ExceptionObject& e )
  {
    std::cerr << e << std::endl;
    exit( -1 );
  }


  return 0;
}

