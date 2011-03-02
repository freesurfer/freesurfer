/**
 * @file  TclProgressDisplayManager.cpp
 * @brief Implements a Tcl based ProgressDisplayManager
 *
 * Uses Tcl functions to implement a progress display. Specifically,
 * uses the NewTask and UpdateTask functions in scuba.tcl which
 * display a dialog with a progress meter.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:39 $
 *    $Revision: 1.10 $
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


#include "string_fixed.h"
#include <sstream>
#include "TclProgressDisplayManager.h"
#include "TclCommandManager.h"

using namespace std;

void
TclProgressDisplayManager::NewTask( string isTitle,
                                    string isText,
                                    bool ibUseMeter,
                                    list<string> ilsButtons ) {

  stringstream ssCommand;
  ssCommand << "NewTask ";
  if ( isTitle != "" ) {
    ssCommand << "-title \"" << isTitle << "\" ";
  }
  if ( isText != "" ) {
    ssCommand << "-text \"" << isTitle << "\" ";
  }
  if ( ibUseMeter ) {
    ssCommand << "-meter true ";
  } else {
    ssCommand << "-meter false ";
  }
  if ( ilsButtons.size() > 0 ) {
    ssCommand << "-buttons {";
    list<string>::iterator tButtons;
    for ( tButtons = ilsButtons.begin(); tButtons != ilsButtons.end();
          ++tButtons ) {
      ssCommand << "\"" << *tButtons << "\" ";
    }
    ssCommand << "}";
  }

  TclCommandManager& manager = TclCommandManager::GetManager();
  manager.SendCommand( ssCommand.str() );

  mlButtons = ilsButtons;

  // Since we're stuck in a c loop at this point, we need to let the
  // Tcl environment handle an event. Otherwise all our windowing
  // stuff will go unnoticed.
  manager.DoTclEvent();

}

void
TclProgressDisplayManager::UpdateTask( string isText,
                                       float iPercent ) {

  stringstream ssCommand;
  ssCommand << "UpdateTask ";
  if ( isText != "" ) {
    ssCommand << "-text \"" << isText << "\" ";
  }
  if ( iPercent != -1 ) {
    ssCommand << "-percent " << iPercent << " ";
  }

  TclCommandManager& manager = TclCommandManager::GetManager();
  manager.SendCommand( ssCommand.str() );

  manager.DoTclEvent();
}

int
TclProgressDisplayManager::CheckTaskForButton() {

  TclCommandManager& manager = TclCommandManager::GetManager();
  string sResult =  manager.SendCommand( "CheckTaskForButtons" );

  // Since we're stuck in a c loop at this point, we need to let the
  // Tcl environment handle an event. Otherwise all our windowing
  // stuff will go unnoticed.
  manager.DoTclEvent();

  // Search the string for the title of each button. If we find it,
  // return its index.
  int nButton = 0;
  list<string>::iterator tButtons;
  for ( tButtons = mlButtons.begin(); tButtons != mlButtons.end();
        ++tButtons ) {

    string::size_type position = sResult.find( *tButtons, 0 );
    if ( position != string::npos ) {
      return nButton;
    }

    nButton++;
  }

  return -1;
}

void
TclProgressDisplayManager::EndTask() {

  TclCommandManager& manager = TclCommandManager::GetManager();
  string sCommand = "EndTask";
  manager.SendCommand( sCommand );
}

