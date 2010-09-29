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

const string inFileDefault = "gcamInput";
#ifdef FS_CUDA
const string outFileDefault = "gpuOutput";
#else
const string outFileDefault = "cpuOutput";
#endif

const float dtDefault = 1;
const float momentumDefault = 0.5;
const float areaIntensityDefault = 0;

string inFilename;
string outFilename;

float dt;
float momentum;
float areaIntensity;

const char* Progname = "gcam_applygradient_test";




// ==========================================================

void ReadCommandLine( int ac, char* av[] ) {

  try {
    bpo::options_description desc("Allowed options");
    desc.add_options()
      ("help", "Produce help message" )
      ("input", bpo::value<string>(&inFilename)->default_value(inFileDefault), "Input gcam filename (.nc will be appended)" )
      ("output", bpo::value<string>(&outFilename)->default_value(outFileDefault), "Output filename (.nc will be appended)" )
      ("dt", bpo::value<float>(&dt)->default_value(dtDefault), "Value of dt")
      ("mom", bpo::value<float>(&momentum)->default_value(momentumDefault), "Valume of momentum")
      ("areaIntensity", bpo::value<float>(&areaIntensity)->default_value(areaIntensityDefault), "Area Intensity" )
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

  cout << "GCAM Apply Gradient Tester" << endl;
  cout << "==========================" << endl << endl;

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
  

  // Set up the GCA_MORPH_PARMS structure
  GCA_MORPH_PARMS parms;
  
  parms.dt = dt;
  parms.momentum = momentum;
  parms.l_area_intensity = areaIntensity;

  // ============================================
  // Perform the calculation
  

#ifdef FS_CUDA
#ifdef GCAMORPH_ON_GPU
  gcamApplyGradientGPU( gcam, &parms );
#else
  cerr << "GCAMORPH_ON_GPU is not defined." << endl;
  cerr << "Test meaningless" << endl;
  exit( EXIT_FAILURE );
#endif
#else
  gcamApplyGradient( gcam, &parms );
#endif

  // ============================================
  // Produce output

  // Fields to check are rx, ry, rz, odx, ody, odz

  myUtils.Write( gcam, outFilename );

  // ====================================
  // Release
  GCAMfree( &gcam );

#ifdef FS_CUDA
  PrintGPUtimers();
#endif

  return( EXIT_SUCCESS );
}
