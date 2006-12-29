/**
 * @file  TclProgressDisplayManager.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 02:09:15 $
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
