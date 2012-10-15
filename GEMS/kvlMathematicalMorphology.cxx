/**
 * @file  kvlMathematicalMorphology.cxx
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
#include "itkOpeningByReconstructionImageFilter.h"
#include "itkMGHImageIOFactory.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImage.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkSubtractImageFilter.h"
#include "itkSimpleFilterWatcher.h"
#include "itkBinaryErodeImageFilter.h"
#include "itkBinaryDilateImageFilter.h"


namespace kvl
{

class ProgressReporter : public itk::SimpleFilterWatcher
{
public:

//
  ProgressReporter( itk::ProcessObject* o, const char *comment="" ) : SimpleFilterWatcher( o, comment )
  {}

protected:

  // Overload some stuff to only show what we need
  virtual void ShowProgress()
  {
    std::cout << "  " << this->GetProcess()->GetProgress() * 100 << "%" << std::endl;
  }

  virtual void StartFilter()
  {
    std::cout << this->GetComment() << ": " << std::endl;
  }

  virtual void EndFilter()
  {
    std::cout << std::endl;
  }

private:

};

};



int main( int argc, char* argv[] )
{
  // Sanity check on input
  if ( argc < 2 )
  {
    std::cerr << argv[0] << " imageFileName [ mode=\"erode\" radius=1 ]" << std::endl;
    std::cerr << std::endl;
    std::cerr << "        with mode is one of \"erode\", \"dilate\", or \"open\"\n\n" << std::endl;
    return -1;
  }

  // Collect command line args
  const std::string  imageFileName = argv[ 1 ];
  std::string  mode = "erode";
  unsigned int  radius = 1;
  if ( argc > 2 )
  {
    mode = argv[ 2 ];
  }
  if ( argc > 3 )
  {
    std::istringstream  radiusStream( argv[ 3 ] );
    radiusStream >> radius;
  }

  std::cout << "imageFileName: " << imageFileName << std::endl;
  std::cout << "radius: " << radius << std::endl;


  // Add support for MGH file format to ITK. An alternative way to add this by default would be
  // to edit ITK's itkImageIOFactory.cxx and explicitly adding it in the code there.
  itk::ObjectFactoryBase::RegisterFactory( itk::MGHImageIOFactory::New() );

  try
  {
    // Set up image reader
    typedef itk::Image< unsigned char, 3 >  ImageType;
    typedef itk::ImageFileReader< ImageType >  ReaderType;
    ReaderType::Pointer  reader = ReaderType::New();
    reader->SetFileName( imageFileName.c_str() );

    // Set up the morphology filter
    typedef itk::BinaryBallStructuringElement< ImageType::PixelType, 3 >  StructuringElementType;
    StructuringElementType  structuringElement;
    structuringElement.SetRadius( radius );
    structuringElement.CreateStructuringElement();

    typedef itk::ImageToImageFilter< ImageType, ImageType >  FilterBaseType;
    FilterBaseType::Pointer  filter = 0;
    if ( !mode.compare( "open" ) )
    {
      std::cout << "Using opening filter" << std::endl;
      typedef itk::OpeningByReconstructionImageFilter< ImageType, ImageType,
              StructuringElementType >  FilterType;
      FilterType::Pointer  myFilter = FilterType::New();
      myFilter->SetKernel( structuringElement );
      filter = myFilter;
    }
    else if ( !mode.compare( "erode" ) )
    {
      std::cout << "Using eroding filter" << std::endl;
      typedef itk::BinaryErodeImageFilter< ImageType, ImageType,
              StructuringElementType >  FilterType;
      FilterType::Pointer  myFilter = FilterType::New();
      myFilter->SetKernel( structuringElement );
      filter = myFilter;
    }
    else if ( !mode.compare( "dilate" ) )
    {
      std::cout << "Using dilating filter" << std::endl;
      typedef itk::BinaryDilateImageFilter< ImageType, ImageType,
              StructuringElementType >  FilterType;
      FilterType::Pointer  myFilter = FilterType::New();
      myFilter->SetKernel( structuringElement );
      filter = myFilter;
    }
    filter->SetInput( reader->GetOutput() );


    // Add an observer
    kvl::ProgressReporter  reporter( filter, mode.c_str() );


    // Set up image writer
    std::ostringstream  outputFileNameStream;
    outputFileNameStream << itksys::SystemTools::GetFilenameWithoutExtension( imageFileName )
                         << "_" << mode << "_radius" << radius << ".mgz";
    const std::string  outputFileName = outputFileNameStream.str();

    typedef itk::ImageFileWriter< ImageType >  WriterType;
    WriterType::Pointer  writer = WriterType::New();
    writer->SetFileName( outputFileName.c_str() );
    writer->SetInput( filter->GetOutput() );

    // Let the beast go!
    writer->Update();

    std::cout << "Wrote output to " << outputFileName << std::endl;
  }
  catch ( itk::ExceptionObject& e )
  {
    std::cerr << e << std::endl;
    exit( -1 );
  }


  return 0;
};

