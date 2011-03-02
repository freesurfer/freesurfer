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


