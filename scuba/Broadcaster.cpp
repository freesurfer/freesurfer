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
 *    $Author: kteich $
 *    $Date: 2007/10/12 22:12:56 $
 *    $Revision: 1.8 $
 *
 * Copyright (C) 2002-2007,
 * The General Hospital Corporation (Boston, MA). 
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 * Bug reports: analysis-bugs@nmr.mgh.harvard.edu
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
