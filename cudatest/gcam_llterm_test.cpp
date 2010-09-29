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

const string inFileDefault = "gcamLLTermInput";
#ifdef FS_CUDA
const string outFileDefault = "gcamLLTermOutputGPU";
#else
const string outFileDefault = "gcamLLTermOutputCPU";
#endif

const string mriFileDefault = "mri.mgz";
const string mriSmoothFileDefault = "mri_smooth.mgz";

const float lambdaDefault = 1;

string inFilename;
string mriFilename;
string mriSmoothFilename;
string outFilename;
float lambda;

const char* Progname = "gcam_llterm_test";




// ==========================================================

void ReadCommandLine( int ac, char* av[] ) {

  try {
    bpo::options_description desc("Allowed options");
    desc.add_options()
      ("help", "Produce help message" )
      ("input", bpo::value<string>(&inFilename)->default_value(inFileDefault), "Input filename (.nc will be appended)" )
      ("mri", bpo::value<string>(&mriFilename)->default_value(mriFileDefault), "MRI filename" )
      ("mrismooth", bpo::value<string>(&mriSmoothFilename)->default_value(mriSmoothFileDefault), "Smoothed MRI filename" )
      ("output", bpo::value<string>(&outFilename)->default_value(outFileDefault), "Output filename (.nc will be appended)" )
      ("lambda", bpo::value<float>(&lambda)->default_value(lambdaDefault), "Value of l_log_likelihood" );
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

  cout << "GCAM LogLikelihood Term Tester" << endl;
  cout << "=============================" << endl << endl;

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

  // Read the MRIs
  MRI* mri = MRIread( mriFilename.c_str() );
  if( !mri ) {
    cerr << "Failed to open " << mriFilename << endl;
    exit( EXIT_FAILURE );
  }

  MRI* mri_smooth = MRIread( mriSmoothFilename.c_str() );
  if( !mri_smooth ) {
     cerr << "Failed to open " << mriSmoothFilename << endl;
     exit( EXIT_FAILURE );
  }

  // ============================================
  // Perform the calculation
  tTotal.Start();
  gcamLogLikelihoodTerm( gcam, mri, mri_smooth, lambda );
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
