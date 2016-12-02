#include "itkImage.h"
#include "itkAddImageFilter.h"
#include "itkMGHImageIOFactory.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"


int main( int argc, char* argv[] )
{
  // Sanity check on input
  if ( argc < 3 )
    {
    std::cerr << argv[0] << " imageFileName1 imageFileName2 [ imageFileName3 ... ]" << std::endl;
    return -1;
    }

  // Add support for MGH file format to ITK. An alternative way to add this by default would be
  // to edit ITK's itkImageIOFactory.cxx and explicitly adding it in the code there.
  itk::ObjectFactoryBase::RegisterFactory( itk::MGHImageIOFactory::New() );

  //
  typedef itk::Image< unsigned char, 3 >  ImageType;

  try
    {
    // Repeatedly read an image, and add it to what we have already
    ImageType::Pointer  result = 0;
    for ( int i = 1; i < argc; i++ )
      {
      // Set up the reader
      typedef itk::ImageFileReader< ImageType >  ReaderType;
      ReaderType::Pointer  reader = ReaderType::New();
      reader->SetFileName( argv[ i ] );

      if ( !result )
        {
        // This is the first image; simply read it
        reader->Update();
        result = reader->GetOutput();
        }
      else
        {
        // We already have something; add that something to the newly read result
        typedef itk::AddImageFilter< ImageType, ImageType, ImageType >  AdderType;
        AdderType::Pointer  adder = AdderType::New();
        adder->SetInput1( result );
        adder->SetInput2( reader->GetOutput() );
        adder->Update();
        result = adder->GetOutput();
        }

      } // End loop over all images


    // Write out the result
    const std::string  resultFileName = "addedProbabilityImage.mgz";
    typedef itk::ImageFileWriter< ImageType >  WriterType;
    WriterType::Pointer  writer = WriterType::New();
    writer->SetInput( result );
    writer->SetFileName( resultFileName.c_str() );
    writer->Update();

    std::cout << "Wrote result to " << resultFileName << std::endl;

    }
  catch( itk::ExceptionObject& e )
    {
    std::cerr << e << std::endl;
    return -1;
    }

  return 0;
};

