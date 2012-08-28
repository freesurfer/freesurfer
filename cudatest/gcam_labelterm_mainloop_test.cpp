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

const string inFileDefault = "gcamLabelTermInput";
#ifdef FS_CUDA
const string outFileDefault = "gcamLabelTermOutputGPU";
#else
const string outFileDefault = "gcamLabelTermOutputCPU";
#endif

const string mriFileDefault = "mri.mgz";
const string mriDistFileDefault = "mri_dist.mgz";
const string gcaFileDefault = "default.gca";

const double lambdaDefault = 1;
const double labelDistDefault = 1;

string inFilename;
string mriFilename;
string mriDistFilename;
string gcaFilename;
string outFilename;
double lambda;
double labelDist;

const char* Progname = "gcam_labelterm__mainloop_test";




// ==========================================================

void ReadCommandLine( int ac, char* av[] ) {

  try {
    bpo::options_description desc("Allowed options");
    desc.add_options()
      ("help", "Produce help message" )
      ("input", bpo::value<string>(&inFilename)->default_value(inFileDefault), "Input filename (.nc will be appended)" )
      ("mri", bpo::value<string>(&mriFilename)->default_value(mriFileDefault), "MRI filename" )
      ("mriDist", bpo::value<string>(&mriDistFilename)->default_value(mriDistFileDefault), "MRI Dist. filename" )
      ("gca", bpo::value<string>(&gcaFilename)->default_value(gcaFileDefault), "GCA filename" )
      ("output", bpo::value<string>(&outFilename)->default_value(outFileDefault), "Output filename (.nc will be appended)" )
      ("lambda", bpo::value<double>(&lambda)->default_value(lambdaDefault), "Value of l_label" )
      ("labelDist", bpo::value<double>(&labelDist)->default_value(labelDistDefault), "Value of label_dist" )
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

  cout << "GCAM Label Term Main Loop Tester" << endl;
  cout << "================================" << endl << endl;

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

  // Read the MRI
  MRI* mri = MRIread( mriFilename.c_str() );
  if( !mri ) {
    cerr << "Failed to open " << mriFilename << endl;
    exit( EXIT_FAILURE );
  }

  // Read mri_dist
  MRI* mri_dist = MRIread( mriDistFilename.c_str() );
  if( !mri_dist ) {
    cerr << "Failed to open " << mriDistFilename << endl;
    exit( EXIT_FAILURE );
  }

  // Read the GCA
  GCA* gca = GCAread( gcaFilename.c_str() );
  if( !gca ) {
    cerr << "Failed to open " << gcaFilename << endl;
    exit( EXIT_FAILURE );
  }

  // Store the gca in the GCAM

  gcam->gca = gca;

  // ============================================
  // Perform the calculation
  tTotal.Start();
  gcamLabelTermMainLoop( gcam, mri, mri_dist, lambda, labelDist ,NULL);
  tTotal.Stop();

  cout << "Computation took " << tTotal << endl;

  // =============================================
  // Fiddle the output
  
  // We need to smuggle out mri_dist
  for( int ix=0; ix<gcam->width; ix++ ) {
    for( int iy=0; iy<gcam->height; iy++ ) {
      for( int iz=0; iz<gcam->depth; iz++ ) {
	gcam->nodes[ix][iy][iz].dx = MRIgetVoxVal(mri_dist,ix,iy,iz,0);
      }
    }
  }

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
