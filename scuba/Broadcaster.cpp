/**
 * @file  Broadcaster.cpp
 * @brief Mixin class that sends messages to Listeners
 *
 * Along with Listener, this implements a broadcaster/listener
 * pattern. Classes can mixin the Broadcaster class and have classes
 * that have mixed in Listener added to their listener list, and then
 * send them all a string message with associated data.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:36 $
 *    $Revision: 1.9 $
 *
 * Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
 *
 * Terms and conditions for use, reproduction, distribution and contribution
 * are found in the 'FreeSurfer Software License Agreement' contained
 * in the file 'LICENSE' found in the FreeSurfer distribution, and here:
 *
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
 *
 * Reporting: freesurfer@nmr.mgh.harvard.edu
 *
 */


#include <iostream>
#include "Broadcaster.h"

using namespace std;

Broadcaster::Broadcaster ( string isLabel ) :
    msLabel(isLabel) {}

Broadcaster::~Broadcaster () {}

void
Broadcaster::AddListener ( Listener& iListener ) {

  mlListeners.push_back( &iListener );
}

void
Broadcaster::RemoveListener ( Listener& iListener ) {

  mlListeners.remove( &iListener );
}

void
Broadcaster::SendBroadcast ( std::string isMessage, void* iData ) {

  // Tell each listener to listen to this message.
  std::list<Listener*>::iterator tListener;
  for ( tListener = mlListeners.begin();
        tListener != mlListeners.end(); ++tListener ) {

    Listener* listener = *tListener;
    listener->ListenToMessage( isMessage, iData );
  }
}

string const&
Broadcaster::GetLabel() const {

  return msLabel;
}
