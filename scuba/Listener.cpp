/**
 * @file  Listener.cpp
 * @brief Mixin class that recieves messages from Broadcasters
 *
 * Along with Broadcaster, this implements a broadcaster/listener
 * pattern. Classes can mixin the Broadcaster class and have classes
 * that have mixed in Listener added to their listener list, and then
 * send them all a string message with associated data.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: kteich $
 *    $Date: 2007/10/12 22:12:56 $
 *    $Revision: 1.6 $
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
Listener::ListenToMessage ( string isMessage, void* iData ) {

  // Call the overrideable function.
  this->DoListenToMessage( isMessage, iData );
}

string const&
Listener::GetLabel() const {

  return msLabel;
}
