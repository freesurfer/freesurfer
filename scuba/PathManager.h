/**
 * @file  PathManager.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 02:09:14 $
 *    $Revision: 1.4 $
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

