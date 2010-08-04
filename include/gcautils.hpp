/**
 * @file  gcautils.hpp
 * @brief C++ GCA utilities
 *
 */
/*
 * Original Authors: Richard Edgar
 * CVS Revision Info:
 *    $Author: rge21 $
 *    $Date: 2010/08/04 19:01:36 $
 *    $Revision: 1.1 $
 *
 * Copyright (C) 2002-2010,
 * The General Hospital Corporation (Boston, MA).
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 *
 */

#ifndef GCA_UTILS_HPP
#define GCA_UTILS_HPP

#include "gca.h"

namespace Freesurfer {

  void GetGCAstats( const GCA* const src );

  void GetGCAnodeStats( const GCA* const src );
  void GetGCApriorStats( const GCA* const src );
}


#endif
