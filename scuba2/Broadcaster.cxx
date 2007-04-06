/**
 * @file  Broadcaster.cxx
 * @brief Mix-in class for message sending
 *
 * Simple mix-in class for use with the Listener class so text
 * messages with a pointer data can be sent to a list of listeners.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: kteich $
 *    $Date: 2007/04/06 22:23:03 $
 *    $Revision: 1.1 $
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
Broadcaster::AddListener ( Listener* const iListener ) {

  mlListeners.push_back( iListener );
}

void
Broadcaster::RemoveListener ( Listener* const iListener ) {

  mlListeners.remove( iListener );
}

void
Broadcaster::SendBroadcast ( std::string const iMessage, 
			     void* const iData ) const {

  //cerr << "Broadcaster " << msLabel << " sending message " << iMessage << endl;

  std::list<Listener*>::const_iterator tListener;
  for ( tListener = mlListeners.begin();
        tListener != mlListeners.end(); ++tListener ) {

    Listener* listener = *tListener;
    // cerr << "\tSending to listener " << listener->GetLabel() << endl;
    listener->ListenToMessage( iMessage, iData );
  }
}
