/**
 * @file  TclScubaKeyCombo.h
 * @brief Tcl specific ScubaKeyCombo
 *
 * This overrides the SetFromString function so we can parse a Tk key
 * string.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: kteich $
 *    $Date: 2007/10/16 16:55:26 $
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


#ifndef TclScubaKeyCombo_h
#define TclScubaKeyCombo_h

#include "ScubaKeyCombo.h"
#include "TclCommandManager.h"

class TclScubaKeyCombo : public ScubaKeyCombo {

  friend class TclScubaKeyComboFactory;
  friend class TclScubaKeyComboStaticTclListener;

public:

  // Overrides constructor to read Tk key strings.
  virtual void SetFromString ( std::string const& isKey );

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

// Our specific factory to create this kind of key type.
class TclScubaKeyComboFactory : public ScubaKeyComboFactory {
public:
  virtual ScubaKeyCombo* MakeKeyCombo() const {
    return new TclScubaKeyCombo();
  }
};



#endif
