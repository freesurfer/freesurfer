/**
 * @file  Listener.cpp
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 02:09:14 $
 *    $Revision: 1.5 $
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
Listener::ListenToMessage ( string iMessage, void* iData ) {

  //  cerr << "Listener " << msLabel << " got message " << iMessage << endl;

  this->DoListenToMessage( iMessage, iData );
}

void
Listener::DoListenToMessage ( string iMessage, void* iData ) {

  this->DoListenToMessage( iMessage, iData );
}


