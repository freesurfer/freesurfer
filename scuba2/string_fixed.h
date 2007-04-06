/**
 * @file  string_fixed.h
 * @brief Wrapper file for including string.h properly
 *
 * Includes a particular version of alloc.h for specific compiler
 * versions. Unfortunately, I can't find the original source of this
 * fix. Please let me know so I can credit it properly.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: kteich $
 *    $Date: 2007/04/06 22:23:04 $
 *    $Revision: 1.1 $
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
