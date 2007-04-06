/**
 * @file  Listener.cxx
 * @brief Mix-in class for message receiving
 *
 * Simple mix-in class for use with the Broadcaster class so text
 * messages with a pointer data can be sent to a list of listeners.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: kteich $
 *    $Date: 2007/04/06 22:23:04 $
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
#include "Listener.h"

using namespace std;

Listener::Listener ( string isLabel ) :
    msLabel( isLabel ) {}

Listener::~Listener () {}

void
Listener::ListenToMessage ( string const iMessage, void* const iData ) {

  //cerr << "Listener " << msLabel << " got message " << iMessage << endl;

  this->DoListenToMessage( iMessage, iData );
}

void
Listener::DoListenToMessage ( string const iMessage, void* const iData ) {

  this->DoListenToMessage( iMessage, iData );
}


