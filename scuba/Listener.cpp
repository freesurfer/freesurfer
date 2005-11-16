#include <iostream>
#include "Listener.h"

using namespace std;

Listener::Listener ( string isLabel ) :
  msLabel( isLabel ) {

}

Listener::~Listener () {
}

void
Listener::ListenToMessage ( string iMessage, void* iData ) {

  //  cerr << "Listener " << msLabel << " got message " << iMessage << endl;
  
  this->DoListenToMessage( iMessage, iData );
}

void
Listener::DoListenToMessage ( string iMessage, void* iData ) {
  
  this->DoListenToMessage( iMessage, iData );
}


