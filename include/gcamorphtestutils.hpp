/**
 * @file  gcamorphtestutils.hpp
 * @brief Utilities to help with testing GCAmorph routines (C++ interface)
 *
 * Reference:
 * "How to Stay Sane while Programming: Collected Wisdom from Broadmoor"
 */
/*
 * Original Author: Richard Edgar
 * CVS Revision Info:
 *    $Author: rge21 $
 *    $Date: 2010/03/11 16:57:23 $
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

#ifndef GCAMORPH_TEST_UTILS_HPP
#define GCAMORPH_TEST_UTILS_HPP


#include <string>
#include <map>

#include <netcdf.h>


#include "gcamorph.h"

// ======================================================================


//! Class to hold utility routines for gcamComputeMetricProperties
class GCAMforCMPutils {
public:

  //! Constructor fills the type map
  GCAMforCMPutils( void );

  //! Writes out relevant portions of a GCAM
  void Write( const GCAM* src, std::string fName ) const;

  //! Reads relevant portions of a GCAM
  void Read( GCAM** dst, string fName ) const;


private:
  //! Map of variable names and types
  std::map<std::string,nc_type> varTypeMap;

  //! Indicies into small arrays defining the dimensions
  enum dimIndices{ iX, iY, iZ };

  //! Number of dimensions we will store
  static const unsigned int nDims = 3;
  //! Number of variables we will store
  static const unsigned int nVars = 8;
};



#endif
