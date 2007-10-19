/**
 * @file  DebugReporter.h
 * @brief Debugging output base class
 *
 * Use this class as a mixin to get basic debug output.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: kteich $
 *    $Date: 2007/10/19 22:31:57 $
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


#ifndef Debug_h
#define Debug_h

#include <iostream>
#include <stdlib.h>

class DebugReporter {

public:
  DebugReporter() {
    mOutputStream = NULL;
    mbDebug = false;

    char* sDebug = getenv( "COMPILE_DEBUG" );
    if ( NULL != sDebug ) {
      mOutputStream = &std::cerr;
      mbDebug = true;
    }
  }

  void DisableOutput () {
    mbDebug = false;
  }

  void SetOutputStream( std::ostream& iStream ) {
    mOutputStream = &iStream;
  }

  void SetOutputStreamToCerr() {
    mOutputStream = &std::cerr;
    mbDebug = true;
  }

#define DebugOutput(x) \
  { \
    mLastReporter = this; \
    if( this->mbDebug ) { \
      if( NULL != mOutputStream ) { \
 (*mOutputStream) << __FILE__ << ":" << __LINE__ << ": " x << std::endl;\
      } else { \
 std::cerr <<__FILE__ << ":" << __LINE__ << ": " x << std::endl; \
      } \
    } \
  }

#define DebugOutputStatic(x) \
  { \
    if( NULL != mLastReporter ) { \
    if( mLastReporter->mbDebug ) { \
      if( NULL != mLastReporter->mOutputStream ) { \
 (*mLastReporter->mOutputStream) << __FILE__ << ":" << __LINE__ << ": " x << std::endl;\
      } else { \
 std::cerr <<__FILE__ << ":" << __LINE__ << ": " x << std::endl; \
      } \
    } \
    } \
  }

  std::ostream* mOutputStream;
  bool mbDebug;

  static DebugReporter* mLastReporter;
};


#endif
