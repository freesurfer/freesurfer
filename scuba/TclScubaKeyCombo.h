#ifndef TclScubaKeyCombo_h
#define TclScubaKeyCombo_h

#include "ScubaKeyCombo.h"
#include "TclCommandManager.h"

class TclScubaKeyCombo : public ScubaKeyCombo {

 public:

  // Overrides constructor to read Tk key strings.
  TclScubaKeyCombo ( std::string );
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

#endif
