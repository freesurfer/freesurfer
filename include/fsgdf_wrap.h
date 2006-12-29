/**
 * @file  fsgdf_wrap.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 02:08:59 $
 *    $Revision: 1.3 $
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


#ifndef fsgdf_wrap_h
#define fsgdf_wrap_h

#include "tcl.h"

#if defined(_WIN32) || defined(__WIN32__)
#       if defined(_MSC_VER)
#               if defined(STATIC_LINKED)
#                       define SWIGEXPORT(a) a
#                       define SWIGIMPORT(a) extern a
#               else
#                       define SWIGEXPORT(a) __declspec(dllexport) a
#                       define SWIGIMPORT(a) extern a
#               endif
#       else
#               if defined(__BORLANDC__)
#                       define SWIGEXPORT(a) a _export
#                       define SWIGIMPORT(a) a _export
#               else
#                       define SWIGEXPORT(a) a
#                       define SWIGIMPORT(a) a
#               endif
#       endif
#else
#       define SWIGEXPORT(a) a
#       define SWIGIMPORT(a) a
#endif

/* -------- TYPES TABLE (END) -------- */

#define SWIG_init    Fsgdf_Init
#define SWIG_name    "fsgdf"
#define SWIG_prefix  ""
#define SWIG_version "0.0"


#ifdef __cplusplus
extern "C"
{
#endif
#ifdef MAC_TCL
#pragma export on
#endif
  SWIGEXPORT(int) SWIG_init(Tcl_Interp *);
#ifdef MAC_TCL
#pragma export off
#endif
#ifdef __cplusplus
}
#endif

#endif
