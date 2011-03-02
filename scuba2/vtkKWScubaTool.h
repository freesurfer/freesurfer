/**
 * @file  vtkKWScubaTool.h
 * @brief Base tool class
 *
 * Base tool class that tools should subclass to implement editing,
 * navigating, etc functionality in reponse to mouse
 * clicks. Subclasses should implement the Do* functions to respond to
 * input events. They can also subclass the AddControls and
 * RemoveControls functions.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:40 $
 *    $Revision: 1.2 $
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


// .NAME vtkKWScubaTool - tool settings for a Scuba window
// .SECTION Description
// Various tool settings and event information. Contains static
// settings such as brush size and target layer, and event information
// like BeginEditEvent, event location. Creates and manages GUI
// widgets for settings.

#ifndef vtkKWScubaTool_h
#define vtkKWScubaTool_h

#include <string>
#include "vtkKWObject.h"
#include "IDTracker.h"

class vtkRenderWindowInteractor;
class vtkKWWidget;
class vtkKWScubaWindow;
class vtkKWScubaView;
class vtkKWScubaLayer;
//BTX
class ScubaInfoItem;
//ETX


class vtkKWScubaTool : public vtkKWObject
      //BTX
      , public IDTracker<vtkKWScubaTool>
      //ETX
{

public:

  static vtkKWScubaTool* New ();
  vtkTypeRevisionMacro( vtkKWScubaTool, vtkKWObject );

  // Description:
  // This should be the label that will be used in the GUI.
  virtual void SetLabel ( const char* isLabel );
  virtual const char* GetLabel ();

  // Description:
  // Populate a UI page with controls for this tool. Depopulate will
  // be called when the tool is no longer displaying its controls on
  // the panel. Subclasses should override AddControls() and
  // RemoveControls().
  void PopulateControlPage ( vtkKWWidget* iPanel );
  void DepopulateControlPage ();

  // Description:
  // This can return true to temporarily suspend picking events during
  // things like click-and-drag redraw events. The navigate tool will
  // do this to retain fast redraws, but the edit tool needs pick
  // events all the time. This also disables mouseover coordinate
  // updating, so use only when necessary.
  virtual bool SuspendPickEvents () { return false; }

  // Description:
  // The view calls these functions, which set internal state
  // information, and then call the internal functions that can be
  // overridden by subclasses to impelement tool behavior.
  void MouseMoveEvent       ( vtkKWScubaWindow* iWindow, vtkKWScubaView* iView,
                              vtkKWScubaLayer* iLayer, float iRAS[3] );
  void LeftMouseDownEvent   ( vtkKWScubaWindow* iWindow, vtkKWScubaView* iView,
                              vtkKWScubaLayer* iLayer, float iRAS[3] );
  void LeftMouseUpEvent     ( vtkKWScubaWindow* iWindow, vtkKWScubaView* iView,
                              vtkKWScubaLayer* iLayer, float iRAS[3] );
  void MiddleMouseDownEvent ( vtkKWScubaWindow* iWindow, vtkKWScubaView* iView,
                              vtkKWScubaLayer* iLayer, float iRAS[3] );
  void MiddleMouseUpEvent   ( vtkKWScubaWindow* iWindow, vtkKWScubaView* iView,
                              vtkKWScubaLayer* iLayer, float iRAS[3] );
  void RightMouseUpEvent    ( vtkKWScubaWindow* iWindow, vtkKWScubaView* iView,
                              vtkKWScubaLayer* iLayer, float iRAS[3] );
  void RightMouseDownEvent  ( vtkKWScubaWindow* iWindow, vtkKWScubaView* iView,
                              vtkKWScubaLayer* iLayer, float iRAS[3] );
  void KeyDownEvent         ( vtkKWScubaWindow* iWindow, vtkKWScubaView* iView,
                              vtkKWScubaLayer* iLayer, float iRAS[3] );
  void EnterEvent           ( vtkKWScubaWindow* iWindow, vtkKWScubaView* iView,
                              vtkKWScubaLayer* iLayer, float iRAS[3] );
  void LeaveEvent           ( vtkKWScubaWindow* iWindow, vtkKWScubaView* iView,
                              vtkKWScubaLayer* iLayer, float iRAS[3] );

protected:

  vtkKWScubaTool ();
  virtual ~vtkKWScubaTool ();

  // Description:
  // Subclasses should override these to add and remove their own
  // controls to the tools panel.
  virtual void AddControls ( vtkKWWidget* iPanel );
  virtual void RemoveControls ();

  // Description:
  // Subclasses should declare these to implement tool behavior. The
  // tool state will be set up, so they can poll themselves. If
  // SuspendPickEvents() returns true, iLayer will be NULL and
  // iRAS will be 0,0,0 for the MouseMoveEvent() callback
  virtual void DoMouseMove  ( vtkKWScubaWindow* iWindow, vtkKWScubaView* iView,
                              vtkKWScubaLayer* iLayer, float iRAS[3] );
  virtual void DoMouseDrag  ( vtkKWScubaWindow* iWindow, vtkKWScubaView* iView,
                              vtkKWScubaLayer* iLayer, float iRAS[3] );
  virtual void DoMouseDown  ( vtkKWScubaWindow* iWindow, vtkKWScubaView* iView,
                              vtkKWScubaLayer* iLayer, float iRAS[3] );
  virtual void DoMouseUp    ( vtkKWScubaWindow* iWindow, vtkKWScubaView* iView,
                              vtkKWScubaLayer* iLayer, float iRAS[3] );
  virtual void DoKeyDown    ( vtkKWScubaWindow* iWindow, vtkKWScubaView* iView,
                              vtkKWScubaLayer* iLayer, float iRAS[3] );
  virtual void DoEnter      ( vtkKWScubaWindow* iWindow, vtkKWScubaView* iView,
                              vtkKWScubaLayer* iLayer, float iRAS[3] );
  virtual void DoLeave      ( vtkKWScubaWindow* iWindow, vtkKWScubaView* iView,
                              vtkKWScubaLayer* iLayer, float iRAS[3] );

  // Description:
  // Called when the tool is clearing its state (e.g. when the mouse
  // leaves the app, and the even is considered done). Subclasses can
  // implement this to clear their own state, or cancel current modal
  // events (e.g. dragging something around).
  virtual void DoClearState ();

  //BTX
  std::string msLabel;
  //ETX

  // Description:
  // Get event information.
  bool IsButtonDown ();
  bool IsButtonUp ();
  int WhichButtonDown ();

  int GetTotalButtonDeltaX ();
  int GetTotalButtonDeltaY ();
  int GetButtonDeltaX ();
  int GetButtonDeltaY ();

  bool IsShiftKeyDown ();
  bool IsControlKeyDown ();
  bool IsAltKeyDown ();

private:

  // Description:
  // Set event information. The tool superclass will do this when the
  // various event callbacks are called.
  void SetStateFromRWI ( vtkRenderWindowInteractor* iRWI );
  void ClearState ();
  void SetWhichButton ( int iButton );
  void SetShiftKeyDown ( bool ibDown );
  void SetControlKeyDown ( bool ibDown );
  void SetAltKeyDown ( bool ibDown );

  int mButton;
  int mTotalDelta[2];
  int mCurrentDelta[2];
  bool mbShiftKeyDown;
  bool mbControlKeyDown;
  bool mbAltKeyDown;

};

#endif
