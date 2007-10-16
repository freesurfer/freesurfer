/**
 * @file  PathManager.h
 * @brief Scuba manager for Path<float>s
 *
 * Managers multiple Path<floats>, relaying pathChanged messages from
 * them to its Listener. Handles Tcl commands to read and write path
 * files.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: kteich $
 *    $Date: 2007/10/16 20:18:30 $
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


#ifndef PathManager_h
#define PathManager_h

#include <vector>

#include "string_fixed.h"
#include "Path.h"
#include "TclCommandManager.h"
#include "Listener.h"
#include "Broadcaster.h"

class PathManager : public TclCommandListener,
      public Broadcaster,  // pathChanged <id>
      public Listener {    // pathChanged <id>

  friend class PathManagerTester;

public:

  // Static function to get our singleton instance.
  static PathManager& GetManager ();

  // While managed, the path will be in the managers list of
  // paths. When unmanaged, it is removed.
  void ManagePath ( Path<float>& iPath );
  void UnmanagePath ( Path<float> const& iPath );

  // Get a list of managed paths.
  std::vector<Path<float>* > const& GetPathList () const;
  std::vector<Path<float>* >& GetPathList ();

  // Respond to Tcl commands.
  virtual TclCommandResult
    DoListenToTclCommand ( char* isCommand, int iArgc, char** iasArgv );

  // On pathChange, passes to listeners.
  virtual void DoListenToMessage ( std::string iMessage, void* iData );

  // Read and write path files.
  void ReadPathFile  ( std::string const& ifnPaths );
  void WritePathFile ( std::string const& ifnPaths );

protected:

  PathManager();

  // Our list of paths.
  std::vector<Path<float>* > mPaths;
};


#endif

