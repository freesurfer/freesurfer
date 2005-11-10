#ifndef Listener_h
#define Listener_h

#include "string_fixed.h"

class Listener {

 public:

  Listener ( std::string isLabel );
  virtual ~Listener ();
  
  void ListenToMessage ( std::string iMessage, void* iData );

  virtual void DoListenToMessage ( std::string iMessage, void* iData );

 protected:
  std::string msLabel;

};

#endif
