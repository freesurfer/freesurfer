/**
 * @file  DebugReporter.h
 * @brief Debugging output base class
 *
 * Use this class as a mixin to get basic debug output.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:36 $
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
