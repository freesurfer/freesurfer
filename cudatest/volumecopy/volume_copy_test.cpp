#include <cstdlib>
#include <cstdio>

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
using namespace std;

#include <boost/program_options.hpp>
namespace bpo = boost::program_options;

#include "mri.h"

#include "chronometer.hpp"

#include "devicemanagement.h"

#include "volcopy.hpp"

// ==========================================================

const string inFilenameDefault = "norm.mgz";
const string outFilenameDefault = "copy.mgz";
const unsigned int kernelSizeDefault = 5;
const unsigned int directionDefault = MRI_WIDTH;



string inFilename;
unsigned int kernelSize;
unsigned int direction;
string outFilename;


const char* Progname = "volume_copy_test";

// ==========================================================

void ReadCommandLine( int ac, char* av[] ) {

  try {
    bpo::options_description desc("Allowed options");
    desc.add_options()
      ("help", "Produce help message" )
      ("input", bpo::value<string>(&inFilename)->default_value(inFilenameDefault), "Input filename" )
      ("output", bpo::value<string>(&outFilename)->default_value(outFilenameDefault), "Output filename" )
      ;
    
    bpo::variables_map vm;
    bpo::store( bpo::parse_command_line( ac, av, desc ), vm );
    bpo::notify( vm );
    
    if( vm.count( "help" ) ) {
      cout << desc << endl;
      exit( EXIT_SUCCESS );
    }
  }
  catch( exception& e ) {
    cerr << "Error: " << e.what() << endl;
    exit( EXIT_FAILURE );
  }
  catch( ... ) {
    cerr << "Unknown exception" << endl;
    exit( EXIT_FAILURE );
  }

} 

// ==========================================================

vector<float> GenerateKernel( const unsigned int nVals ) {
  /*!
    Generates a convolution kernel of the given size
  */

  vector<float> myKernel( nVals );
  float kernelSum = 0;

  // Generate
  for( unsigned int i=0; i<myKernel.size(); i++ ) {
    myKernel.at(i) = i % nVals;
    kernelSum += myKernel.at(i);
  }

  // Normalise
  for( unsigned int i=0; i<myKernel.size(); i++ ) {
    myKernel.at(i) /= kernelSum;
  }

  return( myKernel );
}


// ==========================================================

int main( int argc, char *argv[] ) {


  cout << "VolumeGPU Copy Tester" << endl;
  cout << "=====================" << endl << endl;

  AcquireCUDADevice();

  ReadCommandLine( argc, argv );

  // ======================================
  
  // Read in the input
  cout << "Reading input file: " << inFilename << endl;
  MRI* input = MRIread( inFilename.c_str() );
  if( !input ) {
    cerr << "Failed to open input file: " << inFilename << endl;
    exit( EXIT_FAILURE );
  }

  // Print out some information
  printf( "Sizes are w=%i h=%i d=%i nFrames=%i\n",
	  input->width, input->height, input->depth, input->nframes );
  printf( "Datatype is %i\n\n", input->type );

  // Clone it for output
  MRI* output = NULL;
  output = MRIcopy( input, output );
  


  // ======================================
  // Make sure the output is different from the input

  vector<float> convKernel( GenerateKernel( kernelSizeDefault ) );

  MRIconvolve1d( input, output,
		 &convKernel[0], convKernel.size(),
		 directionDefault,
		 0, 0 );


  // ======================================
  // Do the copy

  VolCopyTest( input, output );

  // ======================================
  // Write out the result

  cout << "Writing output file: " << outFilename << endl;
  MRIwrite( output, outFilename.c_str() );


  return( EXIT_SUCCESS );
}
