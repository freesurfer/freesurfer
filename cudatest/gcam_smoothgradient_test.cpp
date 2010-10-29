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



// ==========================================================

const string inFileDefault = "gcamInput";
#ifdef FS_CUDA
const string outFileDefault = "gcamSmoothGradientOutputGPU";
#else
const string outFileDefault = "gcamSmoothGradientOutputCPU";
#endif

const int nAvgsDefault = 3;

string inFilename;
string outFilename;
int nAvgs;

const char* Progname = "gcam_smoothgradient_test";




// ==========================================================

void ReadCommandLine( int ac, char* av[] ) {

  try {
    bpo::options_description desc("Allowed options");
    desc.add_options()
      ("help", "Produce help message" )
      ("input", bpo::value<string>(&inFilename)->default_value(inFileDefault), "Input filename (.nc will be appended)" )
      ("output", bpo::value<string>(&outFilename)->default_value(outFileDefault), "Output filename (.nc will be appended)" )
      ("navgs", bpo::value<int>(&nAvgs)->default_value(nAvgsDefault), "Size of gaussian smoothing" )
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

  if( nAvgs <= 0 ) {
    std::cerr << "Invalid nAvgs" << std::endl;
    exit( EXIT_FAILURE );
  }
}



// ==========================================================

int main( int argc, char *argv[] ) {

  SciGPU::Utilities::Chronometer tTotal;
  GCAMorphUtils myUtils;

  cout << "GCAM SmoothGradient Tester" << endl;
  cout << "==========================" << endl << endl;

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

  myUtils.Read( &gcam, inFilename );

  // Stop subsequent calls complaining
  gcam->ninputs = 1;


  // ============================================
  // Perform the calculation
  tTotal.Start();
  gcamSmoothGradient( gcam, nAvgs );
  tTotal.Stop();

  cout << "Computation took " << tTotal << endl;


  // =============================================
  // Write the output
  myUtils.Write( gcam, outFilename );

#ifdef FS_CUDA
  PrintGPUtimers();
#endif

  // ====================================
  // Release
  GCAMfree( &gcam );

  exit( EXIT_SUCCESS );
}
