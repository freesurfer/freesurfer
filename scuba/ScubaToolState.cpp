#include <string>
#include "ScubaToolState.h"

using namespace std;

template IDTracker<ScubaToolState>;
int IDTracker<ScubaToolState>::mNextID = 0;
map<int,ScubaToolState*> IDTracker<ScubaToolState>::mIDMap;


ScubaToolState::ScubaToolState() {
  mMode = navigation;

  TclCommandManager& commandMgr = TclCommandManager::GetManager();
  commandMgr.AddCommand( *this, "SetToolMode", 2, "toolID mode",
			 "Sets the current mode of a tool." );
  commandMgr.AddCommand( *this, "GetToolMode", 1, "toolID",
			 "Gets the current mode of a tool." );
}

ScubaToolState::~ScubaToolState() {

}

TclCommandListener::TclCommandResult
ScubaToolState::DoListenToTclCommand ( char* isCommand, 
				       int iArgc, char** iasArgv ) {

  // SetToolMode <toolID> <mode>
  if( 0 == strcmp( isCommand, "SetToolMode" ) ) {
    int toolID = strtol(iasArgv[1], (char**)NULL, 10);
    if( ERANGE == errno ) {
      sResult = "bad tool ID";
      return error;
    }
    
    if( GetID() == toolID ) {

      Mode newMode;
      if( 0 == strcmp( iasArgv[2], "navigation" )) {
	newMode = navigation;
      } else if ( 0 == strcmp( iasArgv[2], "voxelEditing" )) {
	newMode = voxelEditing;
      } else {
	sResult = "bad mode \"" + string(iasArgv[2]) + 
	  "\", should be navigation, voxelEditing.";
	return error;
      }
      SetMode( newMode );
    }
  }

  // GetToolMode <toolID>
  if( 0 == strcmp( isCommand, "GetToolMode" ) ) {
    int toolID = strtol(iasArgv[1], (char**)NULL, 10);
    if( ERANGE == errno ) {
      sResult = "bad tool ID";
      return error;
    }
    
    if( GetID() == toolID ) {

      switch( GetMode() ) {
      case navigation:
	sReturnValues = "navigation";
	break;
      case voxelEditing:
	sReturnValues = "voxelEditing";
	break;
      }
      sReturnFormat = "s";
    }
  }

  return ok;
}
