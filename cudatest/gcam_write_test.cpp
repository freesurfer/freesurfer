#include <cstdlib>
#include <cstdio>

#include <iostream>
#include <iomanip>
#include <string>
using namespace std;

#include <boost/program_options.hpp>
namespace bpo = boost::program_options;

#include "gcamorph.h"
#include "gcamorphtestutils.hpp"

#include "chronometer.hpp"


#ifdef FS_CUDA
#include "devicemanagement.h"
#endif

// ========================================

const string inFileDefault = "gcamInput";
const string outFileDefault = "gcamOutput";

string inFilename; 
string outFilename;

const char* Progname = "gcam_write_test";

// ==========================================================

void ReadCommandLine( int ac, char* av[] ) {

  try {
    bpo::options_description desc("Allowed options");
    desc.add_options()
      ("help", "Produce help message" )
      ("input", bpo::value<string>(&inFilename)->default_value(inFileDefault), "Input filename (.nc will be appended)" )
      ("output", bpo::value<string>(&outFilename)->default_value(outFileDefault), "Output filename (.nc will be appended)" )
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

  GCAM* gcam = NULL;
  GCAMorphUtils myUtils;

  cout << "GCAM Write Tester" << endl;
  cout << "=================" << endl << endl;


  ReadCommandLine( argc, argv );

  // ============================================
  // Read the input file

  myUtils.Read( &gcam, inFilename );

  // =============================================
  // Write the output
  myUtils.Write( gcam, outFilename );

  // ====================================
  // Release
  GCAMfree( &gcam );

#ifdef FS_CUDA
  PrintGPUtimers();
#endif
}
