#include "UndoManager.h"

using namespace std;

UndoManager&
UndoManager::GetManager () {

  static UndoManager* sManager = NULL;
  if( NULL == sManager ) {
    sManager = new UndoManager();

    TclCommandManager& commandMgr = TclCommandManager::GetManager();
    commandMgr.AddCommand( *sManager, "GetUndoTitle", 0, "",
			   "Returns the title for the undo action." );
    commandMgr.AddCommand( *sManager, "UndoOrRedo", 0, "",
			   "Undoes or redoes the current undo list." );
  };

  return *sManager;
}

UndoManager::UndoManager () {
  bUndo = true;
}

void
UndoManager::BeginAction ( std::string isTitle ) {

  // Set our title.
  msTitle = isTitle;

  // Clear our undo list.
  mActions.clear();
}

void
UndoManager::EndAction () {

  
}

string
UndoManager::GetTitle () { 

  if( bUndo ) {
    return "Undo " + msTitle;
  } else {
    return "Redo " + msTitle;
  }
}

void
UndoManager::AddAction ( UndoAction* iAction ) {

  mActions.push_back( iAction );
}

void
UndoManager::Undo () {

  list<UndoAction*>::iterator tActions;
  for( tActions = mActions.begin(); tActions != mActions.end(); ++tActions ) {
    UndoAction* action = *tActions;
    action->Undo();
  }

  bUndo = false;
}

void
UndoManager::Redo () {

  list<UndoAction*>::iterator tActions;
  for( tActions = mActions.begin(); tActions != mActions.end(); ++tActions ) {
    UndoAction* action = *tActions;
    action->Redo();
  }

  bUndo = true;
}

TclCommandManager::TclCommandResult
UndoManager::DoListenToTclCommand ( char* isCommand, 
				    int iArgc, char** iasArgv ) {


  // GetUndoTitle
  if( 0 == strcmp( isCommand, "GetUndoTitle" ) ) {

    sReturnFormat = "s";
    stringstream ssReturnValues;
    ssReturnValues << "\"" << GetTitle() << "\"";
    sReturnValues = ssReturnValues.str();
  }

  // UndoOrRedo
  if( 0 == strcmp( isCommand, "UndoOrRedo" ) ) {
    if( bUndo ) {
      Undo();
    } else {
      Redo();
    }
  }

  return ok;
}

UndoAction::UndoAction () {

}

UndoAction::~UndoAction () {

}

void 
UndoAction::Undo () {

}

void 
UndoAction::Redo () {

}
