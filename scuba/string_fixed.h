/**
 * @file  string_fixed.h
 * @brief Support for older C++ compilers.
 *
 * Supports including the <string> header on old compilers with a
 * weird alloc.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: kteich $
 *    $Date: 2007/10/19 22:31:58 $
 *    $Revision: 1.6 $
 *
 * Copyright (C) 2002-2007,
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


#ifndef string_fixed_h
#define string_fixed_h

#if (__GNUC__ <= 2)
// pick particular version of alloc.h!
#include "/usr/include/g++-3/alloc.h"
#endif

// then
#include <string>


#endif
