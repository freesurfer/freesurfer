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
 *    $Date: 2010/03/09 16:49:35 $
 *    $Revision: 1.1 $
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

  // Construct the filename
  std::string fileName( fName );

  fileName += ".nc";

  std::cout << __FUNCTION__ << ": Writing file " << fileName << std::endl;

  // Reference for the file
  int ncid;

  // Open the file
  NC_SAFE_CALL( nc_create( fileName.c_str(), NC_CLOBBER, &ncid ) );







  // Close the file
  NC_SAFE_CALL( nc_close( ncid ) );
}
