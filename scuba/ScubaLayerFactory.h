#ifndef ScubaLayerFactory_h
#define ScubaLayerFactory_h

#include "string_fixed.h"
#include "DebugReporter.h"
#include "TclCommandManager.h"
#include "Layer.h"

class ScubaLayerFactory : public DebugReporter, public TclCommandListener {
  
  friend class ScubaLayerFactoryTester;

 public:
  // Gets the static reference to this class.
  static ScubaLayerFactory& GetFactory();

  virtual TclCommandResult
    DoListenToTclCommand ( char* iCommand, int iArgc, char** iArgv );

  Layer& MakeLayer ( std::string iType );
  

 protected:
  static bool mbAddedTclCommands;
};

#endif
