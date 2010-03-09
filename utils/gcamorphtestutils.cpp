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
 *    $Date: 2010/03/09 17:45:42 $
 *    $Revision: 1.2 $
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
			    NC_FLOAT, nDims, dimIDs, &varIDs[rz] ) );

  NC_SAFE_CALL( nc_def_var( ncid, "origArea",
			    NC_FLOAT, nDims, dimIDs, &varIDs[origArea] ) );
  NC_SAFE_CALL( nc_def_var( ncid, "area",
			    NC_FLOAT, nDims, dimIDs, &varIDs[area] ) );
  NC_SAFE_CALL( nc_def_var( ncid, "area1",
			    NC_FLOAT, nDims, dimIDs, &varIDs[area2] ) );
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
	
	

	nc_put_vara_double( ncid, varIDs[rx], loc, count, &gcamn.x );
	nc_put_vara_double( ncid, varIDs[ry], loc, count, &gcamn.y );
	nc_put_vara_double( ncid, varIDs[ry], loc, count, &gcamn.z );
      }
    }
  }

  // Close the file
  NC_SAFE_CALL( nc_close( ncid ) );

  std::cout << "complete" << std::endl;
}
