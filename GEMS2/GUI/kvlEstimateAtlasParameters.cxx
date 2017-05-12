#include "kvlAtlasParameterEstimationConsole.h"



int main( int argc, char** argv )
{
  // 
  if ( argc < 2 )
    {
    std::cerr << "Usage: " << argv[ 0 ] << " fileName1 [ fileName2 ... ]" << std::endl;
    exit( -1 );
    }

  // Collect input file names
  std::vector< std::string >  fileNames;
  for ( int i=1; i<argc; i++ )
    {
    fileNames.push_back( std::string( argv[ i ] ) );
    }
  
  //
  kvl::AtlasParameterEstimationConsole  console;
  console.SetLabelImages( fileNames );
  console.Show();

  return 0;
};

