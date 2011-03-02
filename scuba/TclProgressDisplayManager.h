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
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:39 $
 *    $Revision: 1.6 $
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
