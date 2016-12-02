#include "kvlImageThresholdingConsole.h"
#include "itkMGHImageIOFactory.h"


int main( int argc, char* argv[] )
{
  // Sanity check on input
  if ( argc < 2 )
    {
    std::cerr << argv[0] << " imageFileName1 [ imageFileName2 ... ]" << std::endl;
    return -1;
    }

  // Add support for MGH file format to ITK. An alternative way to add this by default would be
  // to edit ITK's itkImageIOFactory.cxx and explicitly adding it in the code there.
  itk::ObjectFactoryBase::RegisterFactory( itk::MGHImageIOFactory::New() );


  // Collect file names
  std::vector< std::string >  fileNames;
  for ( int i = 1; i < argc; i++ )
    {
    fileNames.push_back( argv[ i ] );
    }

  // Set up
  kvl::ImageThresholdingConsole  console;
  console.LoadImages( fileNames );

  // Go!
  console.Show();

  return 0;
};

