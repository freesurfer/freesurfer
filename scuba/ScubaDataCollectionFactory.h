#ifndef ScubaDataCollectionFactory_h
#define ScubaDataCollectionFactory_h

#include <string>
#include "DebugReporter.h"
#include "TclCommandManager.h"
#include "DataCollection.h"

class ScubaDataCollectionFactory : public DebugReporter, public TclCommandListener {
  
  friend class ScubaDataCollectionFactoryTester;

 public:
  // Gets the static reference to this class.
  static ScubaDataCollectionFactory& GetFactory();

  virtual void DoListenToTclCommand ( char* iCommand, int iArgc, char** iArgv );
  DataCollection& MakeDataCollection ( std::string iType );
  

 protected:
  static bool mbAddedTclCommands;
};

#endif
