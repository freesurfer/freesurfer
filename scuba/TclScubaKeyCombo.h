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
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:39 $
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
