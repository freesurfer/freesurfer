#ifndef Listener_h
#define Listener_h

#include <string>

class Listener {

 public:
  
  void ListenToMessage ( std::string iMessage, void* iData );

  virtual void DoListenToMessage ( std::string iMessage, void* iData );

};

#endif
