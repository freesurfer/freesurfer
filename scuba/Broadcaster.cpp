/**
 * @file  Broadcaster.cpp
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 02:09:13 $
 *    $Revision: 1.7 $
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
  for ( tListener = mlListeners.begin();
        tListener != mlListeners.end(); ++tListener ) {

    Listener* listener = *tListener;
    // cerr << "\tSending to listener " << listener->GetLabel() << endl;
    listener->ListenToMessage( iMessage, iData );
  }
}
