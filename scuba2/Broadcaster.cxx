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
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:39 $
 *    $Revision: 1.2 $
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
