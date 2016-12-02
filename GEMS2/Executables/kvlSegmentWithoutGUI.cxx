#include "kvlAtlasMeshSegmentationDriver.h"


int main( int argc, char** argv )
{

  // Sanity check on input
  if ( argc != 2 )
    {
    std::cerr << "Usage: " << argv[ 0 ] << " setUpFileName" << std::endl;
    return( -1 );
    }


  //
  try
    {

    kvl::AtlasMeshSegmentationDriver::Pointer  driver = kvl::AtlasMeshSegmentationDriver::New();
    driver->SetUp( argv[ 1 ] );

    driver->Segment();
    }
  catch( itk::ExceptionObject& e )
    {
    std::cerr << e << std::endl;
    return -1;
    }

  return 0;
};


