#ifndef Debug_h
#define Debug_h

#include <iostream>

class DebugReporter {

 public:
  DebugReporter() {
    mOutputStream = NULL;
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
