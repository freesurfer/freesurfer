/**
 * @file  ScubaROI.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 02:09:14 $
 *    $Revision: 1.6 $
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


#ifndef ScubaROI_h
#define ScubaROI_h

#include "string_fixed.h"
#include "DebugReporter.h"
#include "TclCommandManager.h"
#include "IDTracker.h"
#include "Broadcaster.h"

class ScubaROIStaticTclListener : public DebugReporter,
      public TclCommandListener {

public:
  ScubaROIStaticTclListener ();
  ~ScubaROIStaticTclListener ();

  virtual TclCommandResult
  DoListenToTclCommand ( char* isCommand, int iArgc, char** iasArgv );
};

class ScubaROI : public DebugReporter,
      public IDTracker<ScubaROI>,
      public TclCommandListener,
      public Broadcaster {

  friend class ScubaROITester;

public:

  ScubaROI ();
  virtual ~ScubaROI ();

  virtual TclCommandResult
  DoListenToTclCommand ( char* isCommand, int iArgc, char** iasArgv );

  enum Type { Structure, Free };
  void SetType( Type iType ) {
    mType = iType;
  }
  Type GetType() {
    return mType;
  }

  void SetColorLUT( int iLUTID ) {
    mLUTID = iLUTID;
  }
  int GetColorLUT() {
    return mLUTID;
  }

  void SetStructure( int iStructure ) {
    mStructure = iStructure;
  }
  int GetStructure() {
    return mStructure;
  }

  void SetColor( int iColor[3] );
  void GetColor( int oColor[3] );

  void SetLabel( std::string isLabel ) {
    msLabel = isLabel;
  }
  std::string GetLabel() {
    return msLabel;
  }

  // Returns the color in which this ROI should be drawn; if its type
  // is free, it will return the color, if it's a structure, it will
  // look up its color from an LUT.
  void GetDrawColor( int oColor[3] );

  // Suppresses the roiChanged message. Use when changing a lot of
  // elements in a row that don't need updates in between. Will call
  // ROIChanged() at the end.
  void BeginBatchChanges ();
  void EndBatchChanges ();

protected:

  // For self to call when the ROI had changed. This will broadcast an
  // roiChanged message.
  virtual void ROIChanged ();

  std::string msLabel;

  Type mType;

  int mLUTID;
  int mStructure;
  int mColor[3];

  bool mbSuspendROIChangedMessage;

  static ScubaROIStaticTclListener mStaticListener;
};



#endif
