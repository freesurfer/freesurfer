#ifndef Broadcaster_h
#define Broadcaster_h

#include <list>
#include "string_fixed.h"
#include "Listener.h"

class Broadcaster {

 public:

  Broadcaster ();
  virtual ~Broadcaster ();

  void AddListener ( Listener* iListener );
  void RemoveListener ( Listener* iListener );

  virtual void SendBroadcast ( std::string iMessage, void* iData );

 protected:

  std::list<Listener*> mlListeners;
};


#endif
