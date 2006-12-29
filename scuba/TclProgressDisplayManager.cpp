/**
 * @file  TclProgressDisplayManager.cpp
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 02:09:15 $
 *    $Revision: 1.8 $
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

