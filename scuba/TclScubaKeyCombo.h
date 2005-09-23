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
