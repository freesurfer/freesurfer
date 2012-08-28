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
const string gcaFileDefault = "default.gca";

const double lambdaDefault = 1;
const double labelDistDefault = 1;

string inFilename;
string mriFilename;
string gcaFilename;
string outFilename;
double lambda;
double labelDist;

const char* Progname = "gcam_labelterm_test";




// ==========================================================

void ReadCommandLine( int ac, char* av[] ) {

  try {
    bpo::options_description desc("Allowed options");
    desc.add_options()
      ("help", "Produce help message" )
      ("input", bpo::value<string>(&inFilename)->default_value(inFileDefault), "Input filename (.nc will be appended)" )
      ("mri", bpo::value<string>(&mriFilename)->default_value(mriFileDefault), "MRI filename" )
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

  cout << "GCAM Label Term Tester" << endl;
  cout << "======================" << endl << endl;

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
  gcamLabelTerm( gcam, mri, lambda, labelDist , NULL);
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
