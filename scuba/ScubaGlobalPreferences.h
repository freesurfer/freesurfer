#ifndef ScubaGlobalPreferences_h
#define ScubaGlobalPreferences_h

#include "TclCommandManager.h"

class ScubaGlobalPreferences : public TclCommandListener {

 public:

  // Gets the static reference to this class.
  static ScubaGlobalPreferences& GetPreferences();

  void SavePreferences ();

  virtual TclCommandResult
    DoListenToTclCommand ( char* isCommand, int iArgc, char** iasArgv );

  bool GetViewFlipLeftRightYZ () { return mbViewFlipLeftRightInYZ; }

 protected:

  ScubaGlobalPreferences ();

  void ReadPreferences ();
  void WritePreferences ();

  bool mbViewFlipLeftRightInYZ;
  std::string msInPlaneXKey;
  std::string msInPlaneYKey;
  std::string msInPlaneZKey;
  std::string msCycleKey;
};


#endif
