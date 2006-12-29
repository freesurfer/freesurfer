/**
 * @file  ScubaColorLUT.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 02:09:14 $
 *    $Revision: 1.7 $
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


#ifndef ScubaColorLUT_h
#define ScubaColorLUT_h

#include "string_fixed.h"
extern "C" {
#include "colortab.h"
}
#include "DebugReporter.h"
#include "TclCommandManager.h"
#include "IDTracker.h"

class ScubaColorLUTStaticTclListener : public DebugReporter, public TclCommandListener {

public:
  ScubaColorLUTStaticTclListener ();
  ~ScubaColorLUTStaticTclListener ();

  virtual TclCommandResult
  DoListenToTclCommand ( char* isCommand, int iArgc, char** iasArgv );
};

class ScubaColorLUT : public DebugReporter, public IDTracker<ScubaColorLUT>, public TclCommandListener {

  friend class ScubaColorLUTTester;

public:

  ScubaColorLUT ();
  ~ScubaColorLUT ();

  virtual TclCommandResult
  DoListenToTclCommand ( char* isCommand, int iArgc, char** iasArgv );

  void UseFile ( std::string ifnLUT );

  void GetColorAtIndex ( int iIndex, int oColor[3] );
  void GetIndexOfColor ( int iColor[3], int& oIndex );

  bool IsEntryValid ( int inIndex );

  std::string GetLabelAtIndex ( int iIndex );

  int GetNumberOfEntries () {
    return mEntries.size();
  }

  void SetLabel( std::string isLabel ) {
    msLabel = isLabel;
  }
  std::string GetLabel() {
    return msLabel;
  }

protected:

  static const int cDefaultEntries;

  void ReadFile ();

  std::string mfnLUT;
  int mHighestItemNo;
  struct LUTEntry {
    bool mbValid;
    std::string msLabel;
    int color[3];
    int alpha;
  };
  std::map<int,LUTEntry> mEntries;

  std::string msLabel;

  static ScubaColorLUTStaticTclListener mStaticListener;

};


#endif
