#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"



int main( int argc, char* argv[] )
{

  // Sanity check on input
  if ( argc < 2 )
    {
    std::cerr << "Usage: " << argv[ 0 ] << " erronousManualSegmentation" << std::endl;
    exit( -1 );
    }

  // Collect the input arguments
  const std::string  erronousManualSegmentationFileName( argv[ 1 ] );


  try
    {
    // Read the manual segmentation
    typedef itk::Image< unsigned short, 3 >  ImageType;
    typedef itk::ImageFileReader< ImageType >  ReaderType;
    ReaderType::Pointer  reader = ReaderType::New();
    reader->SetFileName( erronousManualSegmentationFileName.c_str() );
    reader->Update();
    ImageType::Pointer  image = reader->GetOutput();


    // Loop over all voxels, and correct certain values to others
    for ( itk::ImageRegionIterator< ImageType >  it( image, image->GetBufferedRegion() );
          !it.IsAtEnd(); ++it )
      {
      switch ( it.Value() )
        {
        case 550: // left_CA2-3
          it.Value() = 500; // right_CA2-3
          break;
        case 552: // left_CA1  
          it.Value() = 502; // right_CA1
          break;
        case 556: // left_CA4-DG
          it.Value() = 506; // right_CA4-DG
          break;
        }

      } // End loop over all voxels


    // Write out results
    typedef itk::ImageFileWriter< ImageType >  WriterType;
    WriterType::Pointer  writer = WriterType::New();
    writer->SetInput( image );
    writer->SetFileName( "corrected.img" );
    writer->Update();

    }
  catch( itk::ExceptionObject& e )
    {
    std::cerr << e << std::endl;
    }

  return 0;
};
