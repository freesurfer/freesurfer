#include <iostream>
#include "Broadcaster.h"

using namespace std;

Broadcaster::Broadcaster ( string isLabel ) :
  msLabel(isLabel) {

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

  // cerr << "Broadcaster " << msLabel << " sending message " << iMessage << endl;

  std::list<Listener*>::iterator tListener;
  for( tListener = mlListeners.begin(); 
       tListener != mlListeners.end(); ++tListener ) {

    Listener* listener = *tListener;
    // cerr << "\tSending to listener " << listener->GetLabel() << endl;
    listener->ListenToMessage( iMessage, iData );
  }  
}
