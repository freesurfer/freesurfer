/**
 * @file  UndoManager.cpp
 * @brief Manager for a list of undoable (and redoable) actions
 *
 * Use by adding UndoAction subclasses to the list in groups of
 * actions. When told to undo, the Actions will have their Undo()
 * functions called. If an Undo is undone, they will have their Redo()
 * functions called.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:39 $
 *    $Revision: 1.12 $
 *
 * Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
 *
 * Terms and conditions for use, reproduction, distribution and contribution
 * are found in the 'FreeSurfer Software License Agreement' contained
 * in the file 'LICENSE' found in the FreeSurfer distribution, and here:
 *
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
 *
 * Reporting: freesurfer@nmr.mgh.harvard.edu
 *
 */


#ifdef __cplusplus
extern "C"
{
#endif
#include <string.h> // strcmp
#ifdef __cplusplus
}
#endif
#include "UndoManager.h"

using namespace std;

UndoManager&
UndoManager::GetManager () {

  static UndoManager* sManager = NULL;
  if ( NULL == sManager ) {
    sManager = new UndoManager();

    TclCommandManager& commandMgr = TclCommandManager::GetManager();
    commandMgr.AddCommand( *sManager, "GetUndoTitle", 0, "",
                           "Returns the title for the undo action." );
    commandMgr.AddCommand( *sManager, "GetRedoTitle", 0, "",
                           "Returns the title for the redo action." );
    commandMgr.AddCommand( *sManager, "Undo", 0, "",
                           "Undoes the last action." );
    commandMgr.AddCommand( *sManager, "Redo", 0, "",
                           "Redoes last undone action." );
  };

  return *sManager;
}

UndoManager::UndoManager () :
  mcMaxActions( 30 ) {
}

int
UndoManager::BeginAction ( std::string isTitle ) {

  int nNextAction = maCurrentActions.size();

  maCurrentActions[nNextAction] = new UndoableAction();
  maCurrentActions[nNextAction]->msTitle = isTitle;

  return nNextAction;
}

void
UndoManager::AddAction ( int iID, UndoAction* iAction ) {

  if( maCurrentActions.find( iID ) == maCurrentActions.end() ) 
    throw runtime_error( "Tried to add undo action to a non-existant or closed list." );

  maCurrentActions[iID]->mActions.push_back( iAction );
}

void
UndoManager::EndAction ( int iID ) {

  if( maCurrentActions.find( iID ) == maCurrentActions.end() ) 
    throw runtime_error( "Tried to end a non-existant or closed list." );

  mUndoActions.push_back( maCurrentActions[iID] );
  maCurrentActions.erase( iID );
  
  if ( mUndoActions.size() + mRedoActions.size() >
       (unsigned int)mcMaxActions ) {

    // First pop the oldest redo item if we have one, because that's
    // probably the lest neeeded item. If the redo list is empty, do
    // the undo list. The oldest one in each list will be the
    // front. (That's notnecessarily true, but good enough for us.)
    if ( mRedoActions.size() > 0 ) {
      UndoableAction* toDelete = mRedoActions.front();
      mRedoActions.pop_front();
      delete toDelete;
    } else {
      UndoableAction* toDelete = mUndoActions.front();
      mUndoActions.pop_front();
      delete toDelete;
    }
  }
}

string
UndoManager::GetUndoTitle () {

  if ( mUndoActions.size() > 0 ) {
    UndoableAction* undoableAction = mUndoActions.back();
    if ( NULL != undoableAction ) {
      stringstream ss;
      ss << "Undo " << undoableAction->msTitle;
      return ss.str();
    }
  }
  return "No Action to Undo";
}

string
UndoManager::GetRedoTitle () {

  if ( mRedoActions.size() > 0 ) {
    UndoableAction* redoableAction = mRedoActions.back();
    if ( NULL != redoableAction ) {
      stringstream ss;
      ss << "Redo " << redoableAction->msTitle;
      return ss.str();
    }
  }
  return "No Action to Redo";
}

void
UndoManager::Undo () {

  if ( mUndoActions.size() > 0 ) {
    UndoableAction* undoableAction =  mUndoActions.back();
    if ( NULL != undoableAction ) {
      mUndoActions.pop_back();

      // We use a revers iterator here so that the most recently done
      // actions are undone first, and the first done are undone
      // last. This works for an edit where the same voxel is edited
      // multiple times in the same action, as in a brush action, so
      // that the first value is restored last.
      list<UndoAction*>::reverse_iterator tActions;
      for ( tActions = undoableAction->mActions.rbegin();
            tActions != undoableAction->mActions.rend();
            ++tActions ) {

        UndoAction* action = *tActions;
        action->Undo();
      }

      mRedoActions.push_back( undoableAction );
    }
  }
}

void
UndoManager::Redo () {

  if ( mRedoActions.size() > 0 ) {
    UndoableAction* redoableAction = mRedoActions.back();
    if ( NULL != redoableAction ) {
      mRedoActions.pop_back();

      list<UndoAction*>::reverse_iterator tActions;
      for ( tActions = redoableAction->mActions.rbegin();
            tActions != redoableAction->mActions.rend();
            ++tActions ) {
        UndoAction* action = *tActions;
        action->Redo();
      }

      mUndoActions.push_back( redoableAction );
    }
  }
}

void
UndoManager::Clear () {

  list<UndoableAction*>::iterator tAction;
  for ( tAction = mUndoActions.begin();
        tAction != mUndoActions.end();
        ++tAction ) {
    UndoableAction* action = *tAction;
    delete action;
  }

  for ( tAction = mRedoActions.begin();
        tAction != mRedoActions.end();
        ++tAction ) {
    UndoableAction* action = *tAction;
    delete action;
  }

  mUndoActions.clear();
  mRedoActions.clear();
}

TclCommandManager::TclCommandResult
UndoManager::DoListenToTclCommand ( char* isCommand,
                                    int, char** ) {


  // GetUndoTitle
  if ( 0 == strcmp( isCommand, "GetUndoTitle" ) ) {

    sReturnFormat = "s";
    stringstream ssReturnValues;
    ssReturnValues << "\"" << GetUndoTitle() << "\"";
    sReturnValues = ssReturnValues.str();
  }

  // GetRedoTitle
  if ( 0 == strcmp( isCommand, "GetRedoTitle" ) ) {

    sReturnFormat = "s";
    stringstream ssReturnValues;
    ssReturnValues << "\"" << GetRedoTitle() << "\"";
    sReturnValues = ssReturnValues.str();
  }

  // Undo
  if ( 0 == strcmp( isCommand, "Undo" ) ) {
    Undo();
  }

  // Redo
  if ( 0 == strcmp( isCommand, "Redo" ) ) {
    Redo();
  }

  return ok;
}

UndoAction::UndoAction () {}

UndoAction::~UndoAction () {}

void
UndoAction::Undo () {}

void
UndoAction::Redo () {}

UndoableAction::~UndoableAction () {

  list<UndoAction*>::iterator tActions;
  for ( tActions = mActions.begin();
        tActions != mActions.end();
        ++tActions ) {
    UndoAction* action = *tActions;
    delete action;
  }
}
