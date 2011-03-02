/**
 * @file  ScubaColorLUT.h
 * @brief Reads a FreeSurfer LUT file and allows access to values
 *
 * Reads a FressSurfer LUT file (using the utils/colortab.c code) and
 * provides access to values, transforming from colors to indices and
 * getting label strings for indices.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:36 $
 *    $Revision: 1.9 $
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
