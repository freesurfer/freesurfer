#include <iostream>
using namespace std;

#include <boost/program_options.hpp>
namespace bpo = boost::program_options;


#include "chronometer.hpp"

#ifdef FS_CUDA
#include "devicemanagement.h"
#endif


#include "gca.h"
#include "gcalinearnode.hpp"

// ============================================

const string inFileDefault = "input.gca";


string inFilename;

// ============================================

void ReadCommandLine( int ac, char* av[] ) {

  try {
    bpo::options_description desc("Allowed options");
    desc.add_options()
      ("help", "Produce help message" )
      ("input", bpo::value<string>(&inFilename)->default_value(inFileDefault), "Input filename" )
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

  myLinearNode.Exhume( origGCA );
  cout << "Exhumation complete" << endl;
  cout << endl;

  // =============================

  myLinearNode.PrintStats();

  // =============================
  // Free data
  GCAfree( &origGCA );

  return( EXIT_SUCCESS );
}
