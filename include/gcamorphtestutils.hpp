/**
 * @brief Utilities to help with testing GCAmorph routines (C++ interface)
 *
 * Reference:
 * "How to Stay Sane while Programming: Collected Wisdom from Broadmoor"
 */
/*
 * Original Author: Richard Edgar
 *
 * Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
 *
 * Terms and conditions for use, reproduction, distribution and contribution
 * are found in the 'FreeSurfer Software License Agreement' contained
 * in the file 'LICENSE' found in the FreeSurfer distribution, and here:
 *
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
 *
 * Reporting: freesurfer@nmr.mgh.harvard.edu
 *
 */

#ifndef GCAMORPH_TEST_UTILS_HPP
#define GCAMORPH_TEST_UTILS_HPP


#include <string>
#include <map>
#include <iostream>
#include <cstdlib>

#include "netcdf.h"


#include "gcamorph.h"

// ======================================================================

#define NC_SAFE_CALL( call ) do {   \
    int err = call;       \
    if( NC_NOERR != err ) {     \
      std::cerr << __FUNCTION__ \
    << ": NetCDF failure on line " << __LINE__  \
    << " of file " << __FILE__      \
    << std::endl;         \
      std::cerr << "Error code was " << err << std::endl; \
      std::cerr << "Error string was : " << nc_strerror(err)  \
    << std::endl;         \
      abort();                                                  \
    }               \
  } while ( 0 );

// ======================================================================


//! Class to hold utility routines for GCAMorph with one input

class GCAMorphUtils
{
public:

  //! Constructor fills in the type map
  GCAMorphUtils( void );

  //! Writes out a GCAM with one input
  void Write( const GCAM* src, std::string fName ) const;

  //! Reads in a GCAM with one input
  void Read( GCAM** dst, std::string fName ) const;

private:

  //! Indicies into small arrays defining the dimensions
  enum dimIndices { iX, iY, iZ };

  //! Map of variable names and types
  std::map<std::string,nc_type> varTypeMap;

  //! Map of the scalar names and types
  std::map<std::string,nc_type> scalarTypeMap;

  //! Number of dimensions we will store
  static const unsigned int nDims = 3;
  //! Number of variables we will store
  static const unsigned int nVars = 24;
  //! Number of scalars we will store
  static const unsigned int nScalars = 2;

  //! Total number of variables which will be in a NetCDF file
  static const unsigned int totalVars = nVars+nScalars;
};

#endif
