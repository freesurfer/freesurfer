#ifndef ScubaGlobalPreferences_h
#define ScubaGlobalPreferences_h

#include "TclCommandManager.h"

class ScubaGlobalPreferences : public TclCommandListener {

 public:

  enum PrefKey { ShowConsole, ViewFlipLeftRight,
		 KeyInPlaneX, KeyInPlaneY, KeyInPlaneZ, 
		 KeyCycleViewsInFrame };

  // Gets the static reference to this class.
  static ScubaGlobalPreferences& GetPreferences();

  void SavePreferences ();

  virtual TclCommandResult
    DoListenToTclCommand ( char* isCommand, int iArgc, char** iasArgv );

  bool GetPrefAsBool ( PrefKey iKey );

 protected:

  std::string GetStringForKey ( PrefKey iKey );

  ScubaGlobalPreferences ();

  void ReadPreferences ();
  void WritePreferences ();

};


#endif
