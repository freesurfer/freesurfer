/**
 * @file  TclScubaKeyCombo.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 02:09:15 $
 *    $Revision: 1.4 $
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


#ifndef TclScubaKeyCombo_h
#define TclScubaKeyCombo_h

#include "ScubaKeyCombo.h"
#include "TclCommandManager.h"

class TclScubaKeyCombo : public ScubaKeyCombo {

  friend class TclScubaKeyComboFactory;
  friend class TclScubaKeyComboStaticTclListener;

public:

  // Overrides constructor to read Tk key strings.
  virtual void SetFromString ( std::string isKey );

protected:
  TclScubaKeyCombo ();
};

class TclScubaKeyComboStaticTclListener : public TclCommandListener {

public:

  static TclScubaKeyComboStaticTclListener& GetListener ();

  // Registers ConvertTkInputStringToScubaKeyComboString, which takes
  // a string input, makes a TclScubaKeyCombo from it, and returns the
  // ToString() output.
  virtual TclCommandResult
  DoListenToTclCommand ( char* isCommand, int iArgc, char** iasArgv );

protected:

  static bool mbAddedTclCommands;
};

class TclScubaKeyComboFactory : public ScubaKeyComboFactory {
public:
  virtual ScubaKeyCombo* MakeKeyCombo() {
    return new TclScubaKeyCombo();
  }
};



#endif
