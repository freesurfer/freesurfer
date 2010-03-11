/**
 * @file  gcamorphtestutils.cpp
 * @brief Utilities to help with testing GCAmorph routines
 *
 * Reference:
 * "How to Stay Sane while Programming: Collected Wisdom from Broadmoor"
 */
/*
 * Original Author: Richard Edgar
 * CVS Revision Info:
 *    $Author: rge21 $
 *    $Date: 2010/03/11 16:24:45 $
 *    $Revision: 1.11 $
 *
 * Copyright (C) 2002-2008,
 * The General Hospital Corporation (Boston, MA). 
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 * Bug reports: analysis-bugs@nmr.mgh.harvard.edu
 *
 */

#include <string>
#include <vector>
#include <iostream>
#include <map>
using namespace std;

#include <netcdf.h>


#include "gcamorphtestutils.h"

// ======================================================================


#define NC_SAFE_CALL( call ) do {		\
    int err = call; \
    if( NC_NOERR != err ) {			\
      std::cerr << __FUNCTION__ \
		<< ": NetCDF failure on line " << __LINE__	\
		<< " of file " << __FILE__			\
		<< std::endl;					\
      std::cerr << "Error code was " << err << std::endl;	\
      std::cerr << "Error string was : " << nc_strerror(err)	\
		<< std::endl;					\
      exit( EXIT_FAILURE );					\
    }								\
  } while ( 0 );


// ======================================================================

// Class to handle things for gcamComputeMetricProperties test
class GCAMforCMPutils {
public:

  //! Constructor fills the lists
  GCAMforCMPutils( void ) : varTypeMap() {
    
    // Sanity check
    if( this->nVars != 8 ) {
      cerr << __FUNCTION__ << ": Invalid nVars!" << endl;
      exit( EXIT_FAILURE );
    }
    if( this->nDims != 3 ) {
      cerr << __FUNCTION__ << ": Invalid nDims!" << endl;
      exit( EXIT_FAILURE );
    }

    varTypeMap[ "rx" ] = NC_DOUBLE;
    varTypeMap[ "ry" ] = NC_DOUBLE;
    varTypeMap[ "rz" ] = NC_DOUBLE;

    varTypeMap[ "origArea" ] = NC_FLOAT;
    varTypeMap[ "area" ] = NC_FLOAT;
    varTypeMap[ "area1" ] = NC_FLOAT;
    varTypeMap[ "area2" ] = NC_FLOAT;

    varTypeMap[ "invalid" ] = NC_CHAR;

    // And another sanity check
    if( varTypeMap.size() != this->nVars ) {
      cerr << __FUNCTION__ << ": Incorrect entries in varTypeMap" << endl;
      exit( EXIT_FAILURE );
    }
  }



  //! Function to write out relevant portions of a GCAM
  void Write( const GCAM* src, string fName ) const {

    // Construct the filename, using fName passed by value
    fName += ".nc";

    std::cout << __FUNCTION__ << ": Writing file " << fName << " ... ";
    std::cout.flush();

    // Reference for the file
    int ncid;

    // Open the file
    NC_SAFE_CALL( nc_create( fName.c_str(), NC_CLOBBER, &ncid ) );


    // Set up the dimensions
    int dimIDs[nDims];
    NC_SAFE_CALL( nc_def_dim( ncid, "x", src->width, &dimIDs[iX] ) );
    NC_SAFE_CALL( nc_def_dim( ncid, "y", src->height, &dimIDs[iY] ) );
    NC_SAFE_CALL( nc_def_dim( ncid, "z", src->depth, &dimIDs[iZ] ) );

    // Set up the variable IDs
    map<string,int> varIDmap;
    map<string,nc_type>::const_iterator myIt;

    for( myIt = this->varTypeMap.begin();
	 myIt != this->varTypeMap.end();
	 myIt++ ) {
      NC_SAFE_CALL( nc_def_var( ncid,
				myIt->first.c_str(), // Name of the variable
				myIt->second,        // Type of the variable
				nDims, dimIDs,
				&varIDmap[ myIt->first ] ) );
    }

    // Sanity check
    if( varTypeMap.size() != varIDmap.size() ) {
      cerr << __FUNCTION__ << ": Failed to create varIDmap correctly" << endl;
      exit( EXIT_FAILURE );
    }

    // Make the end of the 'definition' region
    NC_SAFE_CALL( nc_enddef( ncid ) );

    // Ugly loop to do the writing element by element
    for( int i=0; i<src->width; i++ ) {
      for( int j=0; j<src->height; j++ ) {
	for( int k=0; k<src->depth; k++ ) {
	  const GCA_MORPH_NODE& gcamn = src->nodes[i][j][k];
	
	  const size_t count[nDims] = { 1, 1, 1};
	  const size_t loc[nDims] = { i, j, k };
	
	
	  // We use 'find' to get an exception if the name doesn't exist
	  NC_SAFE_CALL( nc_put_vara_double( ncid,
					    varIDmap.find( "rx" )->second,
					    loc, count,
					    &gcamn.x ) );
	  NC_SAFE_CALL( nc_put_vara_double( ncid,
					    varIDmap.find( "ry" )->second,
					    loc, count,
					    &gcamn.y ) );
	  NC_SAFE_CALL( nc_put_vara_double( ncid,
					    varIDmap.find( "rz" )->second,
					    loc, count,
					    &gcamn.z ) );

	  NC_SAFE_CALL( nc_put_vara_float( ncid,
					   varIDmap.find( "origArea" )->second,
					   loc, count,
					   &gcamn.orig_area ) );
	  NC_SAFE_CALL( nc_put_vara_float( ncid,
					   varIDmap.find( "area" )->second,
					   loc, count,
					   &gcamn.area ) );
	  NC_SAFE_CALL( nc_put_vara_float( ncid,
					   varIDmap.find( "area1" )->second,
					   loc, count,
					   &gcamn.area1 ) );
	  NC_SAFE_CALL( nc_put_vara_float( ncid,
					   varIDmap.find( "area2" )->second,
					   loc, count,
					   &gcamn.area2 ) );
	  
	  NC_SAFE_CALL( nc_put_vara_text( ncid,
					  varIDmap.find( "invalid" )->second,
					  loc, count,
					  &gcamn.invalid ) );
	}
      }
    }


    // Close the file
    NC_SAFE_CALL( nc_close( ncid ) );

    cout << "complete" << endl;
  }




  //! Reads file containing portions of GCAM used for gcamComputeMetricProperties
  void Read( GCAM** dst, string fName ) const {

    // Make sure input pointer is NULL
    if( *dst != NULL ) {
      cerr << __FUNCTION__
	   << ": dst pointer not NULL!"
	   << endl;
      exit( EXIT_FAILURE );
    }

    // Construct the filename, using fName passed by value
    fName += ".nc";

    std::cout << __FUNCTION__ << ": Reading file " << fName << " ... ";
    std::cout.flush();

    // Reference for the file
    int ncid;

    // Open the file
    NC_SAFE_CALL( nc_open( fName.c_str(), NC_NOWRITE, &ncid ) );

    // Query number of variables and dimensions
    int nDimFile, nVarFile;
    NC_SAFE_CALL( nc_inq_ndims( ncid, &nDimFile ) );
    NC_SAFE_CALL( nc_inq_nvars( ncid, &nVarFile ) );

    if( nDimFile != static_cast<int>(this->nDims) ) {
      cerr << "Invalid number of dimensions " << nDimFile << endl;
      exit( EXIT_FAILURE );
    }

    if( nVarFile != static_cast<int>(this->nVars) ) {
      std::cerr << "Invalid number of variables " << nVarFile << endl;
    }

    int dimIDs[3];
    size_t dimLen[3];

    NC_SAFE_CALL( nc_inq_dimid( ncid, "x", &dimIDs[iX] ) );
    NC_SAFE_CALL( nc_inq_dimid( ncid, "y", &dimIDs[iY] ) );
    NC_SAFE_CALL( nc_inq_dimid( ncid, "z", &dimIDs[iZ] ) );

    NC_SAFE_CALL( nc_inq_dimlen( ncid, dimIDs[iX], &dimLen[iX] ) );
    NC_SAFE_CALL( nc_inq_dimlen( ncid, dimIDs[iY], &dimLen[iY] ) );
    NC_SAFE_CALL( nc_inq_dimlen( ncid, dimIDs[iZ], &dimLen[iZ] ) );

    // Allocate the target
    *dst = GCAMalloc( dimLen[iX], dimLen[iY], dimLen[iZ] );
    if( *dst == NULL ) {
      std::cerr << __FUNCTION__
		<< ": GCAMalloc failed!" << std::endl;
      exit( EXIT_FAILURE );
    }

    // Get hold of the variable IDs
    map<string,int> varIDmap;
    map<string,nc_type>::const_iterator myIt;

    for( myIt = this->varTypeMap.begin();
	 myIt != this->varTypeMap.end();
	 myIt++ ) {
      NC_SAFE_CALL( nc_inq_varid( ncid,
				  myIt->first.c_str(), // Name of the variable
				  &varIDmap[ myIt->first ] ) );
    }

    // Read the data element by element
    for( int i=0; i<(*dst)->width; i++ ) {
      for( int j=0; j<(*dst)->height; j++ ) {
	for( int k=0; k<(*dst)->depth; k++ ) {
	  // Save some type
	  GCA_MORPH_NODE& gcamn = (*dst)->nodes[i][j][k];
	  
	  const size_t count[nDims] = { 1, 1, 1};
	  const size_t loc[nDims] = { i, j, k };
	  
	  
	  NC_SAFE_CALL( nc_get_vara_double( ncid,
					    varIDmap.find( "rx" )->second,
					    loc, count,
					    &gcamn.x ) );
	  NC_SAFE_CALL( nc_get_vara_double( ncid,
					    varIDmap.find( "ry" )->second,
					    loc, count,
					    &gcamn.y ) );
	  NC_SAFE_CALL( nc_get_vara_double( ncid,
					    varIDmap.find( "rz" )->second,
					    loc, count,
					    &gcamn.z ) );

	  NC_SAFE_CALL( nc_get_vara_float( ncid,
					   varIDmap.find( "origArea" )->second,
					   loc, count,
					   &gcamn.orig_area ) );
	  NC_SAFE_CALL( nc_get_vara_float( ncid,
					   varIDmap.find( "area" )->second,
					   loc, count,
					   &gcamn.area ) );
	  NC_SAFE_CALL( nc_get_vara_float( ncid,
					   varIDmap.find( "area1" )->second,
					   loc, count,
					   &gcamn.area1 ) );
	  NC_SAFE_CALL( nc_get_vara_float( ncid,
					   varIDmap.find( "area2" )->second,
					   loc, count,
					   &gcamn.area2 ) );
	  
	  NC_SAFE_CALL( nc_get_vara_text( ncid,
					  varIDmap.find( "invalid" )->second,
					  loc, count,
					  &gcamn.invalid ) );
	}
      }
    }

    
    NC_SAFE_CALL( nc_close( ncid ) );

    cout << "complete" << endl;
  }
  




private:
  //! Map of variable names and types
  map<string,nc_type> varTypeMap;
  enum dimID{ iX, iY, iZ };

  static const unsigned int nDims = 3;
  static const unsigned int nVars = 8;
};



static GCAMforCMPutils myCMPutils;

// ======================================================================

void WriteGCAMforMetricProperties( const GCAM* src, const char* fName ) {
  /*!
    Writes the parts of a GCAM required for gcamComputeMetricProperties
    to a NetCDF file.
    Currently uses the C API, since the C++ one doesn't seem to be available
  */

  myCMPutils.Write( src, fName );
  
}



// ======================================================================

void ReadGCAMforMetricProperties( GCAM** dst, const char* fName ) {
  /*!
    Complement of WriteGCAMforMetricProperties
  */

  myCMPutils.Read( dst, fName );
}
