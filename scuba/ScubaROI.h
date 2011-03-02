/**
 * @file  ScubaROI.h
 * @brief Generic ROI with a an associated structure or color
 *
 * This is a superclass for ROIs. It can be of the type Structure or
 * Free. If Structure, it has a value associated with a ScubaColorLUT
 * and gets a color from that. If it's Free, you can give it an
 * explicit color.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:38 $
 *    $Revision: 1.8 $
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


#ifndef ScubaROI_h
#define ScubaROI_h

#include "string_fixed.h"
#include "DebugReporter.h"
#include "TclCommandManager.h"
#include "IDTracker.h"
#include "Broadcaster.h"

// This is a static listener to respond to static ROI Tcl commands.
class ScubaROIStaticTclListener : public DebugReporter,
      public TclCommandListener {

public:
  ScubaROIStaticTclListener ();
  ~ScubaROIStaticTclListener ();

  // Respond to Tcl commands.
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

  // Respond to Tcl commands.
  virtual TclCommandResult
    DoListenToTclCommand ( char* isCommand, int iArgc, char** iasArgv );
  
  // Our ROI types.
  enum Type { Structure, Free };

  // Get and set the type.
  void SetType( Type iType );
  Type GetType() const;

  // Get and set the ID of the ScubaColorLUT to use if this is a
  // Structure ROI.
  void SetColorLUT( int iLUTID );
  int GetColorLUT() const;

  // Get the structure index.
  void SetStructure( int iStructure );
  int GetStructure() const;

  // Set the free color.
  void SetFreeColor( int const iColor[3] );

  // Human-readable label.
  void SetLabel( std::string const& isLabel );
  std::string const& GetLabel() const;

  // Returns the color in which this ROI should be drawn; if its type
  // is free, it will return the color, if it's a structure, it will
  // look up its color from an LUT.
  void GetDrawColor( int oColor[3] ) const;

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
  int mFreeColor[3];

  bool mbSuspendROIChangedMessage;

  static ScubaROIStaticTclListener mStaticListener;
};



#endif
