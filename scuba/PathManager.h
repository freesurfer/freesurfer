#ifndef PathManager_h
#define PathManager_h


#include "string_fixed.h"
#include <list>
#include "Path.h"
#include "TclCommandManager.h"
#include "Listener.h"
#include "Broadcaster.h"
#include "UndoManager.h"

class PathManager : public TclCommandListener,
		    public Broadcaster,  // pathChanged <id>
                    public Listener {    // pathChanged <id>

  friend class PathManagerTester;

 public:

  static PathManager& GetManager ();

  void ManagePath ( Path<float>& iPath );

  void UnmanagePath ( Path<float>& iPath );

  std::list<Path<float>* >& GetPathList ();

  virtual TclCommandResult
    DoListenToTclCommand ( char* isCommand, int iArgc, char** iasArgv );

  // On pathChange, passes to listeners.
  virtual void DoListenToMessage ( std::string iMessage, void* iData );

  void ReadPathFile  ( std::string ifnPaths );
  void WritePathFile ( std::string ifnPaths );

  void EnableUpdates ();
  void DisableUpdates ();

 protected:
  
  PathManager();

  std::list<Path<float>* > mPaths;
  bool mbSendUpdates;
};


#endif

