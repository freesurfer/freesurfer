/**
 * @file  DebugReporter.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 02:09:13 $
 *    $Revision: 1.5 $
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
