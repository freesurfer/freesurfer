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
    if( this->mbDebug ) { \
      if( NULL != mOutputStream ) { \
	(*mOutputStream) << __FILE__ << ":" << __LINE__ << ": " x << std::endl;\
      } else { \
	std::cerr <<__FILE__ << ":" << __LINE__ << ": " x << std::endl; \
      } \
    } \
  }
  
 protected:

  std::ostream* mOutputStream;
  bool mbDebug;
};


#endif
