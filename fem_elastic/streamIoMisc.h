
#ifndef _h_streamIoMisc_h_
#define _h_streamIoMisc_h_

#include <fstream>

//-------------------
//
// IO Utils
//
//-------------------

template<class T>
void
TWrite(std::ostream& os,
       T value)
{
  T* pt = new T(value);
  os.write( (const char*)pt, sizeof(T) );
  delete pt;
}

template<class T>
T
TRead(std::istream& is)
{
  T* pt = new T;
  is.read( (char*)pt, sizeof(T) );
  T retVal(*pt);
  delete pt;
  return retVal;
}


#endif
