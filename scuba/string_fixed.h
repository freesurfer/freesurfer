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
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:39 $
 *    $Revision: 1.7 $
 *
 * Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
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


#ifndef string_fixed_h
#define string_fixed_h

#if (__GNUC__ <= 2)
// pick particular version of alloc.h!
#include "/usr/include/g++-3/alloc.h"
#endif

// then
#include <string>


#endif
