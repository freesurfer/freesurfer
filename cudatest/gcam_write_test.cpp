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
#include "gcamorphgpu.hpp"
#include "gcamorphcpu.hpp"
#endif

// ========================================

const string inFileDefault = "gcamInput";
const string outFileDefault = "gcamOutput";

string inFilename; 
string outFilename;

#ifdef FS_CUDA
const bool putOnGPUdefault = false;
bool putOnGPU;

const bool linearCPUdefault = false;
bool linearCPU;
#endif

const char* Progname = "gcam_write_test";

// ==========================================================

void ReadCommandLine( int ac, char* av[] ) {

  try {
    bpo::options_description desc("Allowed options");
    desc.add_options()
      ("help", "Produce help message" )
      ("input", bpo::value<string>(&inFilename)->default_value(inFileDefault), "Input filename (.nc will be appended)" )
      ("output", bpo::value<string>(&outFilename)->default_value(outFileDefault), "Output filename (.nc will be appended)" )
#ifdef FS_CUDA
      ("gpu", bpo::value<bool>(&putOnGPU)->default_value(putOnGPUdefault), "Cycle data through GPU" )
      ("linear", bpo::value<bool>(&linearCPU)->default_value(linearCPUdefault), "Cycling data through linear CPU (implies cycle through GPU)" )
#endif
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

  // Optionally cycle through GPU and linear CPU memory
#ifdef FS_CUDA
#ifdef GCAMORPH_ON_GPU
  
  // Coerce putOnGPU if required
  if( linearCPU ) {
    putOnGPU = true;
  }

  if( putOnGPU ) {

    GPU::Classes::GCAmorphGPU myGCAM;

    // Send to the GPU
    myGCAM.SendAll( gcam );

    // Randomise the host allocations
    GPU::Classes::GCAmorphGPU::RandomiseHost();

    if( linearCPU ) {
      Freesurfer::GCAmorphCPU myCPUgcam;

      // Allocate space
      myCPUgcam.AllocateFromTemplate( myGCAM );

      // Get from the GPU
      myCPUgcam.GetFromGPU( myGCAM );

      // Clear the GPU data
      myGCAM.ClearAll();

      // Return to the GPU
      myCPUgcam.PutOnGPU( myGCAM );
    }

    // Retrieve from the GPU
    myGCAM.RecvAll( gcam );
  }
#else
  cerr << "GCAMORPH_ON_GPU is not defined" << endl;
#endif
  cerr << "FS_CUDA is not defined" << endl;
#endif

  // =============================================
  // Write the output
  myUtils.Write( gcam, outFilename );

  // ====================================
  // Release
  GCAMfree( &gcam );

#ifdef FS_CUDA
  PrintGPUtimers();
#ifdef GCAMORPH_ON_GPU
  GPU::Classes::GCAmorphGPU::ReleaseHost();
#endif
#endif

  return( EXIT_SUCCESS );
}
