/**
 * @file  ProgressDisplayManager.h
 * @brief Generic progress update GUI functions
 *
 * This is a virtual class meant to be reimplemented in a specific
 * widget framework. See TclProgressDisplayManager and
 * QtProgressDisplayManager. Clients can use this class to start a
 * task and check for user input during that class, i.e. to cancel it.
 */
/*
 * Original Author:Kevin Teich
 * CVS Revision Info:
 *    $Author: kteich $
 *    $Date: 2007/10/22 04:39:27 $
 *    $Revision: 1.3 $
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


#ifndef ProgressDisplayManager_h
#define ProgressDisplayManager_h

#include "string_fixed.h"
#include <list>

class ProgressDisplayManager {

public:

  // Use this to set a specific implementation of the class, i.e.:
  // ProgressDisplayManager::SetManager( new TclProgressDisplayManager );
  static void SetManager ( ProgressDisplayManager* iManager );

  // Get the manager.
  static ProgressDisplayManager& GetManager();

  ProgressDisplayManager() {}
  virtual ~ProgressDisplayManager() {}

  // Start a new task with a title and text explanation. ibUseMeter
  // specifies to use an update meter; use it if you know how long
  // your task will take and can update it at regular
  // intervals. ilsButtons is a list of text buttons.  Override this
  // and put up a dialog box or whatever with buttons.
  virtual void NewTask ( std::string isTitle,
                         std::string isText,
                         bool ibUseMeter,
                         std::list<std::string> ilsButtons ) = 0;

  // Use this if you are using a progress bar. Pass a 0-100 percent
  // and a text string.  In subclass, update the progress bar.
  virtual void UpdateTask ( std::string isText,
                            float iPercent ) = 0;

  // If the user pressed on of your buttons, this will return the
  // index of it. Will return -1 if no buttons were pressed.  In
  // subclass, you should be checking for button inputs. Return the
  // index of a button that was pressed, according to the order passed
  // in from NewTask.
  virtual int CheckTaskForButton () = 0;

  // Please call this function when the task is ended.  Subclasses
  // should use this to tear down dialog boxes.
  virtual void EndTask () = 0;

private:

  static ProgressDisplayManager* sManager;
};


#endif
