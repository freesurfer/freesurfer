#include <iostream>
#include <string>
using namespace std;

#include <boost/program_options.hpp>
namespace bpo = boost::program_options;


#include "chronometer.hpp"

#ifdef FS_CUDA
#include "devicemanagement.h"
#endif


#include "gca.h"
#include "gcalinearnode.hpp"
#include "gcalinearprior.hpp"

// ============================================

const string inFileDefault = "input.gca";
const string outFileDefault = "output.gca";


string inFilename;
string outFilename;

// ============================================

void ReadCommandLine( int ac, char* av[] ) {

  try {
    bpo::options_description desc("Allowed options");
    desc.add_options()
      ("help", "Produce help message" )
      ("input", bpo::value<string>(&inFilename)->default_value(inFileDefault), "Input filename" )
      ("output", bpo::value<string>(&outFilename)->default_value(outFileDefault), "Input filename" )
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



// ============================================

const char* Progname = "gca_read_test";

int main( int argc, char *argv[] ) {
  
  cout << "GCA Read Test" << endl;
  cout << "=============" << endl;
  cout << endl;

  
  ReadCommandLine( argc, argv );

  // =============================
  // Read in the GCA

  GCA* origGCA = GCAread( inFilename.c_str() );
  cout << "Read in GCA from file" << endl;
  cout << endl;


  // =============================
  // Pack GCA into linear format
  Freesurfer::GCAlinearNode myLinearNode;
  Freesurfer::GCAlinearPrior myLinearPrior;

  myLinearNode.Exhume( origGCA );
  myLinearPrior.Exhume( origGCA );
  cout << "Exhumation complete" << endl;

  // =============================
  // Send back to the GCA
  myLinearNode.Inhume( origGCA );
  myLinearPrior.Inhume( origGCA );
  cout << "Inhumation complete" << endl;

  // =============================
  // Write the output
  GCAwrite( origGCA, outFilename.c_str() );
  cout << "Written output file" << endl;
  

  // =============================
  // Free data
  GCAfree( &origGCA );

  // =============================
  // Print the timers
  myLinearNode.PrintStats();
  myLinearPrior.PrintStats();

  return( EXIT_SUCCESS );
}
