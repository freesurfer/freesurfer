#ifndef ScubaToolState_h
#define ScubaToolState_h

#include "IDTracker.h"
#include "TclCommandManager.h"

class ScubaToolState : public TclCommandListener, public IDTracker<ScubaToolState> {

 public:
  ScubaToolState();
  virtual ~ScubaToolState();

  enum Mode { navigation, voxelEditing };

  void SetMode ( Mode iMode ) { mMode = iMode; }
  Mode GetMode () { return mMode; }

  virtual TclCommandResult
    DoListenToTclCommand ( char* iCommand, int iArgc, char** iasArgv );

 protected:

  Mode mMode;

};


#endif
