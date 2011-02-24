#include <cstdlib>
#include <cstdio>

#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
using namespace std;

#ifdef TRAP_NAN
#include <fenv.h>
#endif

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


const int spacingDefault = 1;
int spacing;

string inFilename;
string outFilename;


const char* Progname = "gcam_removesingularities_test";




// ==========================================================

void ReadCommandLine( int ac, char* av[] ) {

  try {
    bpo::options_description desc("Allowed options");
    desc.add_options()
      ("help", "Produce help message" )
      ("input", bpo::value<string>(&inFilename)->default_value(inFileDefault), "Input gcam filename (.nc will be appended)" )
      ("output", bpo::value<string>(&outFilename)->default_value(outFileDefault), "Output filename (.nc will be appended)" )
      ("spacing", bpo::value<int>(&spacing)->default_value(spacingDefault), "Spacing value in gcam" )
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

#ifdef TRAP_NAN
  feenableexcept( FE_INVALID );
#endif

  cout << "GCAM Remove Singularities Tester" << endl;
  cout << "================================" << endl << endl;

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
  
  // ============================================
  // Perform the calculation
  
  gcam->spacing = spacing;

  tTotal.Start();
#ifdef FS_CUDA
#ifdef GCAMORPH_ON_GPU
  GCAMremoveSingularitiesGPU( gcam );
#else
  cerr << "GCAMORPH_ON_GPU is not defined." << endl;
  cerr << "Test meaningless" << endl;
  exit( EXIT_FAILURE );
#endif
#else
  MRI *mri_warp = NULL;
  mri_warp = GCAMwriteWarpToMRI(gcam, mri_warp);
  printf( "mri_warp : %i %i %i\n",
          mri_warp->width, mri_warp->height, mri_warp->depth );
  printf( "gcam : %i %i %i\n",
          gcam->width, gcam->height, gcam->depth );
  GCAMremoveSingularitiesAndReadWarpFromMRI( gcam, mri_warp );
  MRIfree( &mri_warp);
#endif
  tTotal.Stop();
  cout << "Computation took " << tTotal << endl;

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
