#ifndef UndoManager_h
#define UndoManager_h

#include <string>
#include <list>
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

class UndoManager : public TclCommandListener {

  friend class UndoManagerTester;

 public:

  static UndoManager& GetManager ();

  // Call to begin and end a new undoable action. In between call
  // AddAction with subclassed UndoActions to implement the
  // action. The title should be something that works with "Undo" or
  // "Redo" at the beginning. It will show up in the Edit menu bar. So
  // use something like "Selection", so it comes out as "Undo
  // Selection" and "Redo Selection."
  void BeginAction ( std::string isTitle );
  void EndAction ();

  std::string GetTitle ();

  // Adds an action to the undo list.
  void AddAction ( UndoAction* iAction );

  // Calls Undo or Redo on all the actions in our undo list.
  void Undo ();
  void Redo ();

  virtual TclCommandResult
    DoListenToTclCommand ( char* isCommand, int iArgc, char** iasArgv );

 protected:
  
  UndoManager();

  std::string msTitle;

  // List of undo actions.
  std::list<UndoAction*> mActions;

  // If true, we're ready to undo. If false, we're ready to redo.
  bool bUndo;
};


#endif 
