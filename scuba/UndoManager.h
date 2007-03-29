/**
 * @file  UndoManager.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: kteich $
 *    $Date: 2007/03/29 21:36:35 $
 *    $Revision: 1.5 $
 *
 * Copyright (C) 2002-2007,
 * The General Hospital Corporation (Boston, MA). 
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 * Bug reports: analysis-bugs@nmr.mgh.harvard.edu
 *
 */


#ifndef UndoManager_h
#define UndoManager_h

#include "string_fixed.h"
#include <list>
#include <map>
#include "TclCommandManager.h"


// Subclass this class to implement undoable actions.
class UndoAction {

public:

  UndoAction ();
  virtual ~UndoAction ();

  // These will be called when the action is to be undone or
  // redone. For example, for an UndoAction that undoes the selection
  // of a voxel, Undo() should unselect that voxel, and Redo() should
  // reselect the voxel. The subclass should have all the necessary
  // information to perform this action.
  virtual void Undo ();
  virtual void Redo ();
};

class UndoableAction {

public:

  UndoableAction() {}
  virtual ~UndoableAction ();

  std::list<UndoAction*> mActions;
  std::string msTitle;

};

class UndoManager : public TclCommandListener {

  friend class UndoManagerTester;

public:

  static UndoManager& GetManager ();
  
  // Call to begin and end a new undoable action. In between call
  // AddAction with subclassed UndoActions to implement the
  // action. The title should be something that works with "Undo" or
  // "Redo" at the beginning. It will show up in the Edit menu bar. So
  // use something like "Selection", so it comes out as "Undo
  // Selection" and "Redo Selection." Beginning an action returns an
  // ID which should be passed to AddAction and EndAction to specify
  // the action.
  int BeginAction ( std::string isTitle );
  void EndAction ( int iID );

  // Adds an action to the Undo list.
  void AddAction ( int iID, UndoAction* iAction );

  // Get titles for our undo and redo actions.
  std::string GetUndoTitle ();
  std::string GetRedoTitle ();

  // Calls Undo or Redo on all the actions in our undo list. Calling
  // Undo takes an action off the undo list and puts it on the Redo
  // list, and vice versa.
  void Undo ();
  void Redo ();

  virtual TclCommandResult
  DoListenToTclCommand ( char* isCommand, int iArgc, char** iasArgv );

  // Clear undo and redo lists.
  void Clear();

protected:

  UndoManager();

  // Actions we're currently building.
  std::map<int,UndoableAction*> maCurrentActions;

  // List of undo and redo actions.
  std::list<UndoableAction*> mUndoActions;
  std::list<UndoableAction*> mRedoActions;

  int mcMaxActions;
};


#endif
