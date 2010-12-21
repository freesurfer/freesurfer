#include <cstdlib>
#include <cstdio>

#include <string>
#include <iomanip>
#include <iostream>
using namespace std;


#include <boost/program_options.hpp>
namespace bpo = boost::program_options;

#include "mri.h"

#include "chronometer.hpp"

#ifdef FS_CUDA
#include "devicemanagement.h"
#endif

// ========================================

const string inFilenameDefault = "mri.mgz";

#ifdef FS_CUDA
const string outFilenameStemDefault = "gpu";
#else
const string outFilenameStemDefault = "cpu";
#endif

const int markDefault = 1;
const int labelStartDefault = 0;
const int labelStopDefault = 40;
const int sixConnectDefault = 0;


string inFilename;
string outFilenameStem;
int mark;
int labelStart;
int labelStop;
int sixConnect;


const char* Progname = "mriMarkLabelBorderVoxelTest";


// ==========================================================

void ReadCommandLine( int ac, char* av[] ) {

  try {
    bpo::options_description desc("Allowed options");
    desc.add_options()
      ("help", "Produce help message" )
      ("input", bpo::value<string>(&inFilename)->default_value(inFilenameDefault), "Input filename" )
      ("outstem", bpo::value<string>(&outFilenameStem)->default_value(outFilenameStemDefault), "Stem for output filenames" )
      ("mark", bpo::value<int>(&mark)->default_value(markDefault), "Value for marks" )
      ("start", bpo::value<int>(&labelStart)->default_value(labelStartDefault), "First label to examine" )
      ("stop", bpo::value<int>(&labelStop)->default_value(labelStopDefault), "Last label to examine" )
      ("sixConnect", bpo::value<int>(&sixConnect)->default_value(sixConnectDefault), "Six connect test" )
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

int main( int argc, char *argv[] ) {

  SciGPU::Utilities::Chronometer tTotal;

  cout << "MRImarkLabelBorderVoxel Tester" << endl;
  cout << "==============================" << endl << endl;

  
#ifdef FS_CUDA
  AcquireCUDADevice();
#else
  cout << "CPU Version" << endl;
#endif
  cout << endl;

  ReadCommandLine( argc, argv );
  
  // ----------------------------------------------

  // Read in the input
  cout << "Reading input file: " << inFilename << endl;
  MRI* input = MRIread( inFilename.c_str() );
  if( !input ) {
    cerr << "Failed to open input file: " << inFilename << endl;
    exit( EXIT_FAILURE );
  }

  // ----------------------------------------------
  for( int i=labelStart; i<=labelStop; i++ ) {

    MRI* output = NULL;

    tTotal.Start();
    output = MRImarkLabelBorderVoxels( input, output,
				       i, mark, sixConnect );
    tTotal.Stop();

    stringstream ofname;
    ofname << outFilenameStem
	   << setw(3) << setfill('0') << i
	   << ".mgz";
    cout << "Writing file: " << ofname.str() << endl;

    MRIwrite( output, ofname.str().c_str() );
    
  }

  cout << "Compute time: " << tTotal << endl;

  // ----------------------------------------------

#ifdef FS_CUDA
  PrintGPUtimers();
#endif

  MRIfree( &input );

  return( EXIT_SUCCESS );
}
