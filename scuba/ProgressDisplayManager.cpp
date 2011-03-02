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
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:36 $
 *    $Revision: 1.4 $
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

