/**
 * @file  ProgressDisplayManager.cpp
 * @brief Generic progress update GUI functions
 *
 * This is a virtual class meant to be reimplemented in a specific
 * widget framework. See TclProgressDisplayManager and
 * QtProgressDisplayManager. Clients can use this class to start a
 * task and check for user input during that class, i.e. to cancel it.
 */
/*
 * Original Author:Kevin Teich
 * CVS Revision Info:
 *    $Author: kteich $
 *    $Date: 2007/10/22 04:39:27 $
 *    $Revision: 1.3 $
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


#include <stdexcept>
#include "ProgressDisplayManager.h"

ProgressDisplayManager* ProgressDisplayManager::sManager = NULL;

void
ProgressDisplayManager::SetManager ( ProgressDisplayManager* iManager ) {

  sManager = iManager;
}


ProgressDisplayManager&
ProgressDisplayManager::GetManager() {

  if ( NULL == sManager ) {
    throw std::runtime_error( "No manager set yet." );
  }

  return *sManager;
}

