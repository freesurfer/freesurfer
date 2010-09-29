#include <cstdlib>
#include <cstdio>

#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
using namespace std;

#include <boost/program_options.hpp>
namespace bpo = boost::program_options;

#include "gcamorph.h"
#include "gcamorphtestutils.hpp"

#include "chronometer.hpp"


#ifdef FS_CUDA
#include "devicemanagement.h"
#endif



// ==========================================================

const string gcamFileDefault = "gcamInput";
#ifdef FS_CUDA
const string outFileDefault = "energyLabel.gpu";
#else
const string outFileDefault = "energyLabel.cpu";
#endif


string gcamFilename;
string outFilename;


const char* Progname = "gcam_labelenergy_test";




// ==========================================================

void ReadCommandLine( int ac, char* av[] ) {

  try {
    bpo::options_description desc("Allowed options");
    desc.add_options()
      ("help", "Produce help message" )
      ("gcam", bpo::value<string>(&gcamFilename)->default_value(gcamFileDefault), "Input gcam filename (.nc will be appended)" )
      ("output", bpo::value<string>(&outFilename)->default_value(outFileDefault), "Output filename" )
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
  GCAMorphUtils myUtils;

  cout << "GCAM Label Energy Tester" << endl;
  cout << "========================" << endl << endl;


#ifdef FS_CUDA
#ifndef GCAMORPH_ON_GPU
  cerr << "GCAMORPH_ON_GPU is not defined." << endl;
  cerr << "Test meaningless" << endl;
  exit( EXIT_FAILURE );
#endif
#endif


#ifdef FS_CUDA
  AcquireCUDADevice();
#else
  cout << "CPU Version" << endl;
#endif
  cout << endl;

  ReadCommandLine( argc, argv );

  // ============================================

  // Read the input file

  GCAM* gcam = NULL;
  myUtils.Read( &gcam, gcamFilename );
  
  // ============================================
  // Perform the calculation
  tTotal.Start();
  double labelenergy = gcamLabelEnergy( gcam, NULL, 0 );
  tTotal.Stop();

  cout << "Computation took " << tTotal << endl;


  // ============================================
  // Produce output
  cout << "Label Energy = " << labelenergy << endl;

  cout << "Writing file " << outFilename << endl;

  ofstream oFile;
  oFile.open( outFilename.c_str(), ios::binary | ios::out | ios::trunc );

  oFile.write( reinterpret_cast<const char*>(&labelenergy),
	       sizeof(labelenergy) );
  oFile.close();

  // ====================================
  // Release
  GCAMfree( &gcam );

#ifdef FS_CUDA
  PrintGPUtimers();
#endif
  
  exit( EXIT_SUCCESS );
}
