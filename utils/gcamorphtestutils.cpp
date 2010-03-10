/**
 * @file  gcamorphtestutils.c
 * @brief Utilities to help with testing GCAmorph routines
 *
 * Reference:
 * "How to Stay Sane while Programming: Collected Wisdom from Broadmoor"
 */
/*
 * Original Author: Richard Edgar
 * CVS Revision Info:
 *    $Author: rge21 $
 *    $Date: 2010/03/10 18:52:52 $
 *    $Revision: 1.6 $
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
#include <iostream>

#include <netcdf.h>


#include "gcamorphtestutils.h"

const unsigned int nDims = 3;
enum dimID{ iX, iY, iZ };

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

void WriteGCAMforMetricProperties( const GCAM* src, const char* fName ) {
  /*!
    Writes the parts of a GCAM required for gcamComputeMetricProperties
    to a NetCDF file.
    Currently uses the C API, since the C++ one doesn't seem to be available
  */

  enum myVars{ rx, ry, rz, origArea, invalid, area, area1, area2 };

  if( src->ninputs != 1 ) {
    std::cerr << __FUNCTION__ << ": Must have ninputs=1!" << std::endl;
    exit( EXIT_FAILURE );
  }

  // Construct the filename
  std::string fileName( fName );

  fileName += ".nc";

  std::cout << __FUNCTION__ << ": Writing file " << fileName << " ... ";
  std::cout.flush();

  // Reference for the file
  int ncid;

  // Open the file
  NC_SAFE_CALL( nc_create( fileName.c_str(), NC_CLOBBER, &ncid ) );


  // Set up the dimension references
  int dimIDs[nDims];
  NC_SAFE_CALL( nc_def_dim( ncid, "x", src->width, &dimIDs[iX] ) );
  NC_SAFE_CALL( nc_def_dim( ncid, "y", src->height, &dimIDs[iY] ) );
  NC_SAFE_CALL( nc_def_dim( ncid, "z", src->depth, &dimIDs[iZ] ) );

  

  // Set up the variables
  int varIDs[8];
  NC_SAFE_CALL( nc_def_var( ncid, "r_x",
			    NC_DOUBLE, nDims, dimIDs, &varIDs[rx] ) );
  NC_SAFE_CALL( nc_def_var( ncid, "r_y",
			    NC_DOUBLE, nDims, dimIDs, &varIDs[ry] ) );
  NC_SAFE_CALL( nc_def_var( ncid, "r_z",
			    NC_DOUBLE, nDims, dimIDs, &varIDs[rz] ) );

  NC_SAFE_CALL( nc_def_var( ncid, "origArea",
			    NC_FLOAT, nDims, dimIDs, &varIDs[origArea] ) );
  NC_SAFE_CALL( nc_def_var( ncid, "area",
			    NC_FLOAT, nDims, dimIDs, &varIDs[area] ) );
  NC_SAFE_CALL( nc_def_var( ncid, "area1",
			    NC_FLOAT, nDims, dimIDs, &varIDs[area1] ) );
  NC_SAFE_CALL( nc_def_var( ncid, "area2",
			    NC_FLOAT, nDims, dimIDs, &varIDs[area2] ) );


  NC_SAFE_CALL( nc_def_var( ncid, "invalid",
			    NC_CHAR, nDims, dimIDs, &varIDs[invalid] ) );


  // Make the end of the 'definition' region
  NC_SAFE_CALL( nc_enddef( ncid ) );


  // Write the data (in a very ugly way)
  for( int i=0; i<src->width; i++ ) {
    for( int j=0; j<src->height; j++ ) {
      for( int k=0; k<src->depth; k++ ) {
	const GCA_MORPH_NODE& gcamn = src->nodes[i][j][k];
	
	const size_t count[nDims] = { 1, 1, 1};
	const size_t loc[nDims] = { i, j, k };
	
	

	NC_SAFE_CALL( nc_put_vara_double( ncid, varIDs[rx],
					  loc, count, &gcamn.x ) );
	NC_SAFE_CALL( nc_put_vara_double( ncid, varIDs[ry],
					  loc, count, &gcamn.y ) );
	NC_SAFE_CALL( nc_put_vara_double( ncid, varIDs[ry],
					  loc, count, &gcamn.z ) );
	
	NC_SAFE_CALL( nc_put_vara_float( ncid, varIDs[origArea],
					 loc, count, &gcamn.orig_area ) );
	NC_SAFE_CALL( nc_put_vara_float( ncid, varIDs[area],
					 loc, count, &gcamn.area ) );
	NC_SAFE_CALL( nc_put_vara_float( ncid, varIDs[area1],
					 loc, count, &gcamn.area1 ) );
	NC_SAFE_CALL( nc_put_vara_float( ncid, varIDs[area2],
					 loc, count, &gcamn.area2 ) );

	NC_SAFE_CALL( nc_put_vara_text( ncid, varIDs[invalid],
					loc, count, &gcamn.invalid ) );
      }
    }
  }

  // Close the file
  NC_SAFE_CALL( nc_close( ncid ) );

  std::cout << "complete" << std::endl;
}



// ======================================================================

void ReadGCAMforMetricProperties( GCAM** dst, const char* fName ) {
  /*!
    Complement of WriteGCAMforMetricProperties
  */

  if( *dst != NULL ) {
    std::cerr << __FUNCTION__
	      << ": dst pointer not NULL!"
	      << std::endl;
    exit( EXIT_FAILURE );
  }
 

  // Construct the filename
  std::string fileName( fName );

  fileName += ".nc";

  std::cout << __FUNCTION__ << ": Reading file " << fileName << std::endl;

  // Reference for the file
  int ncid;

  // Open the file
  NC_SAFE_CALL( nc_open( fileName.c_str(), NC_NOWRITE, &ncid ) );

  // Query number of variables and dimensions
  int nDim, nVar;
  NC_SAFE_CALL( nc_inq_ndims( ncid, &nDim ) );
  NC_SAFE_CALL( nc_inq_nvars( ncid, &nVar ) );

  if( nDim != 3 ) {
    std::cerr << "Invalid number of dimensions " << nDim << std::endl;
    exit( EXIT_FAILURE );
  }

  if( nVar != 8 ) {
    std::cerr << "Invalid number of variables " << nVar << std::endl;
  }

  // Fetch the dimensions
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
  enum myVars{ rx, ry, rz, origArea, invalid, area, area1, area2 };
  int varIDs[8];

  NC_SAFE_CALL( nc_inq_varid( ncid, "r_x", &varIDs[rx] ) );
  NC_SAFE_CALL( nc_inq_varid( ncid, "r_y", &varIDs[ry] ) );
  NC_SAFE_CALL( nc_inq_varid( ncid, "r_z", &varIDs[rz] ) );

  NC_SAFE_CALL( nc_inq_varid( ncid, "origArea", &varIDs[origArea] ) );
  NC_SAFE_CALL( nc_inq_varid( ncid, "area", &varIDs[area] ) );
  NC_SAFE_CALL( nc_inq_varid( ncid, "area1", &varIDs[area1] ) );
  NC_SAFE_CALL( nc_inq_varid( ncid, "area2", &varIDs[area2] ) );

  NC_SAFE_CALL( nc_inq_varid( ncid, "invalid", &varIDs[invalid] ) );

  // Read the data element by element
  for( int i=0; i<(*dst)->width; i++ ) {
    for( int j=0; j<(*dst)->height; j++ ) {
      for( int k=0; k<(*dst)->depth; k++ ) {
	// Save some type
	GCA_MORPH_NODE& gcamn = (*dst)->nodes[i][j][k];

	const size_t count[nDims] = { 1, 1, 1};
	const size_t loc[nDims] = { i, j, k };

	NC_SAFE_CALL( nc_get_vara_double( ncid, varIDs[rx],
					  loc, count, &gcamn.x ) );
	NC_SAFE_CALL( nc_get_vara_double( ncid, varIDs[ry],
					  loc, count, &gcamn.y ) );
	NC_SAFE_CALL( nc_get_vara_double( ncid, varIDs[rz],
					  loc, count, &gcamn.z ) );

	NC_SAFE_CALL( nc_get_vara_float( ncid, varIDs[origArea],
					 loc, count, &gcamn.orig_area ) );
	NC_SAFE_CALL( nc_get_vara_float( ncid, varIDs[area],
					 loc, count, &gcamn.area ) );
	NC_SAFE_CALL( nc_get_vara_float( ncid, varIDs[area1],
					 loc, count, &gcamn.area1 ) );
	NC_SAFE_CALL( nc_get_vara_float( ncid, varIDs[area2],
					 loc, count, &gcamn.area2 ) );

	NC_SAFE_CALL( nc_get_vara_text( ncid, varIDs[invalid],
					loc, count, &gcamn.invalid ) );
	
      }
    }
  }


  NC_SAFE_CALL( nc_close( ncid ) );

  std::cout << __FUNCTION__ << ": Read complete" << std::endl;
}
