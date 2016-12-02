#include "kvlAtlasParameterEstimationConsole.h"



int main( int argc, char** argv )
{
  // 
  if ( argc < 4 )
    {
    std::cerr << "Usage: " << argv[ 0 ] << " useGaussians ignoreLastImage fileName1 [ fileName2 ... ]" << std::endl;
    exit( -1 );
    }

  // Collect useGaussians
  std::istringstream  useGaussiansStream( argv[ 1 ] );
  bool  useGaussians;
  useGaussiansStream >> useGaussians;

  // Collect ignoreLastImage
  std::istringstream  ignoreLastImageStream( argv[ 2 ] );
  bool  ignoreLastImage;
  ignoreLastImageStream >> ignoreLastImage;


  // Collect input file names
  std::vector< std::string >  fileNames;
  for ( int i=3; i<argc; i++ )
    {
    fileNames.push_back( std::string( argv[ i ] ) );
    }
  
  //
  kvl::AtlasParameterEstimationConsole  console;
  console.SetLabelImages( fileNames, useGaussians, ignoreLastImage );
  console.Show();

  return 0;
};

