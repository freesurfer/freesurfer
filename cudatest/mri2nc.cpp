/*! \file
  This is a program to concatenate multiple MRI scans into a single
  NetCDF file.
  The scans all have to have the same dimensions, but they do not
  have to be of the same type.
  For now, only frame 0 is taken of each scan.
*/
#include <cstdlib>

#include <string>
#include <iostream>
#include <vector>
#include <map>
using namespace std;

#include <netcdf.h>

#include <boost/program_options.hpp>
namespace bpo = boost::program_options;


#include "mri.h"


// ======================================================================

#define NC_SAFE_CALL( call ) do {               \
    int err = call;                             \
    if( NC_NOERR != err ) {                     \
      std::cerr << __FUNCTION__					\
                << ": NetCDF failure on line " << __LINE__      \
                << " of file " << __FILE__                      \
                << std::endl;                                   \
      std::cerr << "Error code was " << err << std::endl;       \
      std::cerr << "Error string was : " << nc_strerror(err)    \
                << std::endl;                                   \
      exit( EXIT_FAILURE );                                     \
    }                                                           \
  } while ( 0 );

const unsigned int nDims = 3;

// This ordering makes the x index run fastest
enum dimIndices{ iZ, iY, iX };



// ==========================================================


const string outFileDefault = "mri.nc";


vector<string> mriFiles;
string outFile;

const char* Progname = "mri2nc";


// ==========================================================

void ReadCommandLine( int ac, char* av[] ) {

  try {
    
    bpo::options_description generic("Generic options");
    generic.add_options()
      ("help", "Produce help message" )
      ("version,v", "Show source version" )
      ;

    bpo::options_description output( "Output options" );
    output.add_options()
      ("output", bpo::value<string>(&outFile)->default_value(outFileDefault), "Name of the output file" )
      ;
    

    bpo::options_description files("Hidden input files");
    files.add_options()
      ("inputFiles", bpo::value< vector<string> >(&mriFiles), "List of MRI files" )
      ;

    bpo::positional_options_description p;
    p.add( "inputFiles", -1 );

    bpo::options_description cmdLine;
    cmdLine.add(files).add(generic).add(output);

    
    bpo::options_description visible( "Allowed options" );
    visible.add(generic).add(output);
    

    bpo::variables_map vm;
    bpo::store( bpo::command_line_parser(ac, av).
		options(cmdLine).positional(p).run(), vm );
    bpo::notify( vm );

    if( vm.count( "help" ) ) {
      cout << "Usage :" << endl;
      cout << "  mri2nc [--output=<file.nc>] <mri1> <mri2> <mri3> ..." << endl;
      cout << visible << endl;
      exit( EXIT_SUCCESS );
    }

    if( vm.count( "version" ) ) {
      cout << __FILE__ << endl
	   << "$Author: rge21 $\n"
	   << "$Date: 2010/04/29 16:08:36 $\n"
	   << "$Revision: 1.1 $\n" 
	   << endl;
      exit( EXIT_SUCCESS );
    }

    if( !vm.count( "inputFiles" ) ) {
      cerr << "No input files specified!" <<endl;
      exit( EXIT_FAILURE );
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

void LoadMRIfiles( const vector<string> fileNames,
		   map<string,MRI*>& mriData ) {
  vector<string>::const_iterator file;

  for( file=fileNames.begin(); file!=fileNames.end(); file++ ) {
    mriData[ *file ] = MRIread( file->c_str() );
    if( !mriData[ *file ] ) {
      cerr << "Failed to open " << *file << endl;
      exit( EXIT_FAILURE );
    }
  }
}


// ==========================================================

void VerifyDims( const map<string,MRI*>& mriData ) {

  map<string,MRI*>::const_iterator mriIt=mriData.begin();

  // Get the first values
  const int width = mriIt->second->width;
  const int height = mriIt->second->height;
  const int depth = mriIt->second->depth;
  
  mriIt++;
  
  while( mriIt != mriData.end() ) {
    if( (width  != (mriIt->second->width)) ||
	(height != (mriIt->second->height)) ||
	(depth  != (mriIt->second->depth)) ) {
      cerr << "Dimension mismatch in "
	   << mriIt->first << endl;
      exit( EXIT_FAILURE );
    }

    mriIt++;
  }

  cout << "Dimensions are :"
       << "( " << width
       << ", " << depth
       << ", " << height << " )" << endl;
    
  // If this exits, check was OK
}

// ==========================================================

template<typename T>
void ExhumeRow( const MRI* src,
		T* h_row,
		const unsigned int iy,
		const unsigned int iz,
		const unsigned int iFrame ) {
  cerr << __FUNCTION__
       << ": Unrecognised data type "
       << src->type << endl;
  exit( EXIT_FAILURE );
  h_row[0] = iy * iz * iFrame;
}


template<>
void ExhumeRow<unsigned char>( const MRI* src,
			       unsigned char* h_row,
			       const unsigned int iy,
			       const unsigned int iz,
			       const unsigned int iFrame ) {
  memcpy( h_row,
	  &MRIseq_vox( src, 0, iy, iz, iFrame ),
	  src->width*sizeof(unsigned char) );
}

template<>
void ExhumeRow<float>( const MRI* src,
		       float* h_row,
		       const unsigned int iy,
		       const unsigned int iz,
		       const unsigned int iFrame ) {
  memcpy( h_row,
	  &MRIFseq_vox( src, 0, iy, iz, iFrame ),
	  src->width*sizeof(float) );
}


template<typename T>
void ExhumeMRI( const MRI* src,
		T* data ) {
  // Reimplementation of ExhumeFrame method of MRIframeGPU
  const unsigned int iFrame = 0;

  unsigned int currLoc = 0;

  for( int iz=0; iz<src->depth; iz++ ) {
    for( int iy=0; iy<src->height; iy++ ) {
      
      ExhumeRow( src, &(data[currLoc]), iy, iz, iFrame );

      // Add the offset for the next row
      currLoc += src->width;
    }
  }
}


// ==========================================================

string TrimVarName( const string fullName ) {

  string res;

  size_t foundDir = fullName.find_last_of( "/" );
  size_t foundExt = fullName.find_last_of( "." );

  res = fullName.substr( foundDir+1, foundExt );

  cout << fullName << " became " << res << endl;
  return( res );
}

// ==========================================================

void WriteNetCDF( const map<string,MRI*>& mriData,
		  const string fName ) {
  
  int ncid;
  map<string,MRI*>::const_iterator mriIt=mriData.begin();
  map<string,int> varIDmap;

  // Open the file
  NC_SAFE_CALL( nc_create( fName.c_str(), NC_CLOBBER, &ncid ) );

  // Set up the dimensions
  int dimIDs[nDims];
  
  // Have already checked that dimensions match
  NC_SAFE_CALL( nc_def_dim( ncid, "x", mriIt->second->width, &dimIDs[iX] ) );
  NC_SAFE_CALL( nc_def_dim( ncid, "y", mriIt->second->height, &dimIDs[iY] ) );
  NC_SAFE_CALL( nc_def_dim( ncid, "z", mriIt->second->depth, &dimIDs[iZ] ) );

  for( mriIt = mriData.begin();
       mriIt != mriData.end();
       mriIt++ ) {

    // Get the 'normalised' variable name
    string varName = TrimVarName( mriIt->first );

    switch( mriIt->second->type ) {
    case MRI_UCHAR:
      NC_SAFE_CALL( nc_def_var( ncid,
				varName.c_str(),
				NC_BYTE,
				nDims, dimIDs,
				&varIDmap[ mriIt->first ] ) );
      break;

    case MRI_FLOAT:
      NC_SAFE_CALL( nc_def_var( ncid,
				varName.c_str(),
				NC_FLOAT,
				nDims, dimIDs,
				&varIDmap[ mriIt->first ] ) );
      break;

    default:
      cerr << __FUNCTION__
	   << ": Unrecognised data type in MRI "
	   << mriIt->first
	   << " on line " << __LINE__ << endl;
      exit( EXIT_FAILURE );
    }
  }

  // Mark the end of the 'definition' region
  NC_SAFE_CALL( nc_enddef( ncid ) );

  // Calculate the number of elements
  mriIt = mriData.begin();
  const size_t nElems = mriIt->second->width
    * mriIt->second->height
    * mriIt->second->depth;

  for( mriIt=mriData.begin();
       mriIt!=mriData.end();
       mriIt++ ) {

    vector<unsigned char> dataChar( nElems );
    vector<float> dataFloat( nElems );

    switch( mriIt->second->type ) {

    case MRI_UCHAR:
      ExhumeMRI( mriIt->second, &(dataChar[0]) );
      NC_SAFE_CALL( nc_put_var_uchar( ncid,
				      varIDmap.find( mriIt->first )->second,
				      &dataChar[0] ) );
      break;

    case MRI_FLOAT:
      ExhumeMRI( mriIt->second, &(dataFloat[0]) );
      NC_SAFE_CALL( nc_put_var_float( ncid,
				      varIDmap.find( mriIt->first )->second,
				      &dataFloat[0] ) );
      break;

    default:
      cerr << __FUNCTION__
	   << ": Unrecognised data type in MRI "
	   << mriIt->first
	   << " on line " << __LINE__ << endl;
      exit( EXIT_FAILURE );
    }

  }


  // Close the file
  NC_SAFE_CALL( nc_close( ncid ) );
}



// ===========================================================

int main( int argc, char *argv[] ) {

  // Get command line arguments
  ReadCommandLine( argc, argv );

  // Load the files
  map<string,MRI*> mriData;

  LoadMRIfiles( mriFiles, mriData );

  VerifyDims( mriData );

  WriteNetCDF( mriData, outFile );

  return( EXIT_SUCCESS );
}
