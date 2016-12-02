#include "kvlImageViewingConsole.h"
#include "itkMGHImageIOFactory.h"


int main( int argc, char* argv[] )
{
  // Sanity check on input
  if ( argc < 2 )
    {
    std::cerr << argv[0] << " imageFileName [ overlayImageFileName overlayLookupTable meshCollectionFileName meshNumber ]" << std::endl;
    return -1;
    }

  // Add support for MGH file format to ITK. An alternative way to add this by default would be
  // to edit ITK's itkImageIOFactory.cxx and explicitly adding it in the code there.
  itk::ObjectFactoryBase::RegisterFactory( itk::MGHImageIOFactory::New() );


  //
  kvl::ImageViewingConsole  console;
  console.LoadImage( argv[ 1 ] );

  if ( argc > 2 )
    {
    console.LoadOverlayImage( argv[ 2 ] );
    }

  if ( argc > 3 )
    {
    console.LoadOverlayLookupTable( argv[ 3 ] );
    }

  if ( argc > 4 )
    {
    int  meshNumber = -1;
    if ( argc > 5 )
      {
      std::istringstream  meshNumberStream( argv[ 5 ] );
      meshNumberStream >> meshNumber;
      }
    
    console.LoadMesh( argv[ 4 ], meshNumber );
    }

  console.Show();

  return 0;
};

