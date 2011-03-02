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
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:36 $
 *    $Revision: 1.7 $
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
