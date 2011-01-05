#include <cstdlib>
#include <cstdio>

#include <string>
#include <iomanip>
#include <iostream>
using namespace std;


#include <boost/program_options.hpp>
namespace bpo = boost::program_options;

#include "mri.h"

#include "chronometer.hpp"

#ifdef FS_CUDA
#include "devicemanagement.h"
#endif

// ========================================

const string mriFilenameDefault = "mri.mgz";
const string mrivalsFilenameDefault = "mri_vals.mgz";
const string mixFilenameDefault = "mix";
const string nbrFilenameDefault = "nbr";

#ifdef FS_CUDA
const string outFilenameStemDefault = "gpu";
#else
const string outFilenameStemDefault = "cpu";
#endif

const int labelStartDefault = 0;
const int labelStopDefault = 40;


string mriFilename;
string mrivalsFilename;
string mixFilename;
string nbrFilename;
string outFilenameStem;
int labelStart;
int labelStop;


const char* Progname = "mriVoxInLabelPartialVolumeTest";


// ==========================================================

void ReadCommandLine( int ac, char* av[] ) {

  try {
    bpo::options_description desc("Allowed options");
    desc.add_options()
      ("help", "Produce help message" )
      ("mri", bpo::value<string>(&mriFilename)->default_value(mriFilenameDefault), "Input mri filename" )
      ("mrivals", bpo::value<string>(&mrivalsFilename)->default_value(mrivalsFilenameDefault), "Input mri_vals filename" )
      ("outstem", bpo::value<string>(&outFilenameStem)->default_value(outFilenameStemDefault), "Stem for output filenames" )
      ("mix", bpo::value<string>(&mixFilename)->default_value(mixFilenameDefault), "Middle of 'mix' MRI filename" )
      ("nbr", bpo::value<string>(&nbrFilename)->default_value(nbrFilenameDefault), "Middle of 'nbr' MRI filename" )
      ("start", bpo::value<int>(&labelStart)->default_value(labelStartDefault), "First label to examine" )
      ("stop", bpo::value<int>(&labelStop)->default_value(labelStopDefault), "Last label to examine" )
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


// ------------------------

void WriteMRI( const MRI* mri,
	       const string stem,
	       const string id,
	       const int label ) {

  stringstream ofname;

  ofname << stem
	 << id
	 << setw(3) << setfill('0') << label
	 << ".mgz";
  
  MRIwrite( const_cast<MRI*>(mri), ofname.str().c_str() );
}


// ==========================================================

int main( int argc, char *argv[] ) {

  SciGPU::Utilities::Chronometer tTotal;

  cout << "mriVoxInLabelPartialVolume Tester" << endl;
  cout << "=================================" << endl << endl;

  
#ifdef FS_CUDA
  AcquireCUDADevice();
#else
  cout << "CPU Version" << endl;
#endif
  cout << endl;

  ReadCommandLine( argc, argv );
  
  // ----------------------------------------------

  // Read in the input
  cout << "Reading input file: " << mriFilename << endl;
  MRI* mri = MRIread( mriFilename.c_str() );
  if( !mri ) {
    cerr << "Failed to open input file: " << mriFilename << endl;
    exit( EXIT_FAILURE );
  }

  cout << "Reading input file: " << mrivalsFilename << endl;
  MRI* mri_vals = MRIread( mrivalsFilename.c_str() );
  if( !mri_vals ) {
    cerr << "Failed to open input file: " << mrivalsFilename << endl;
    exit( EXIT_FAILURE );
  }

  // ----------------------------------------------
  for( int i=labelStart; i<=labelStop; i++ ) {

    MRI* mix = NULL;
    MRI* nbr = NULL;

    mix = MRIalloc( mri->width, mri->height, mri->depth, MRI_FLOAT );
    nbr = MRIalloc( mri->width, mri->height, mri->depth, MRI_UCHAR );

    float vol;

    tTotal.Start();
    vol =MRIvoxelsInLabelWithPartialVolumeEffects( mri,
						   mri_vals,
						   i,
						   mix,
						   nbr );
    tTotal.Stop();


    cout << "Volume was " << vol << endl;

    // Create some output files
    WriteMRI( nbr, outFilenameStem, nbrFilename, i );
    WriteMRI( mix, outFilenameStem, mixFilename, i );
    

    MRIfree( &mix );
    MRIfree( &nbr );
  }

  cout << "Compute time: " << tTotal << endl;

  // ----------------------------------------------

#ifdef FS_CUDA
  PrintGPUtimers();
#endif

  MRIfree( &mri );
  MRIfree( &mri_vals );

  return( EXIT_SUCCESS );
}
