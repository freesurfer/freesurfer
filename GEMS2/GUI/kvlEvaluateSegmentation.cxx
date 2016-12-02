#include "kvlSegmentationEvaluationConsole.h"


int main( int argc, char* argv[] )
{
  // Sanity check on input
  if ( argc < 4 )
    {
    std::cerr << argv[0] << " imageFileName overlayImageFileName1 overlayImageFileName2 [ overlayLookupTable ]" << std::endl;
    return -1;
    }
  
  //
  kvl::SegmentationEvaluationConsole  console;
  console.LoadImage( argv[ 1 ] );
  console.LoadOverlayImage1( argv[ 2 ] );
  console.LoadOverlayImage2( argv[ 3 ] );

  if ( argc > 4 )
    {
    console.LoadOverlayLookupTable( argv[ 4 ] );
    }

  console.Show();

  return 0;
};

