#include <iostream>
#include "Broadcaster.h"

using namespace std;

Broadcaster::Broadcaster () {

}

Broadcaster::~Broadcaster () {

}

void
Broadcaster::AddListener ( Listener* iListener ) {

  mlListeners.push_back( iListener );
}

void
Broadcaster::RemoveListener ( Listener* iListener ) {

  mlListeners.remove( iListener );
}

void
Broadcaster::SendBroadcast ( std::string iMessage, void* iData ) {

  std::list<Listener*>::iterator tListener;
  for( tListener = mlListeners.begin(); 
       tListener != mlListeners.end(); ++tListener ) {

    Listener* listener = *tListener;
    listener->ListenToMessage( iMessage, iData );
  }  
}
