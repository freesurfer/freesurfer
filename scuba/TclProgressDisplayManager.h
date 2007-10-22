/**
 * @file  TclProgressDisplayManager.h
 * @brief Implements a Tcl based ProgressDisplayManager
 *
 * Uses Tcl functions to implement a progress display. Specifically,
 * uses the NewTask and UpdateTask functions in scuba.tcl which
 * display a dialog with a progress meter.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: kteich $
 *    $Date: 2007/10/22 04:39:30 $
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


#ifndef TclProgressDisplayManager_h
#define TclProgressDisplayManager_h

#include "ProgressDisplayManager.h"
#include "string_fixed.h"
#include <list>

class TclProgressDisplayManager : public ProgressDisplayManager {

public:

  TclProgressDisplayManager() {}
  ~TclProgressDisplayManager() {}

  void NewTask ( std::string isTitle,
                 std::string isText,
                 bool ibUseMeter,
                 std::list<std::string> ilsButtons );

  void UpdateTask ( std::string isText,
                    float iPercent );

  int CheckTaskForButton ();

  void EndTask ();

protected:

  std::list<std::string> mlButtons;

};


#endif
