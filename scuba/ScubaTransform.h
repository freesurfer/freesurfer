/**
 * @file  ScubaTransform.h
 * @brief A Transform44 object wrapped for Scuba
 *
 * This is a Transform44 that behaves like a Scuba object, with
 * IDTrackeing, Tcl scripting, and Broadcasting of its changes. It
 * also manages the special case of when a transform is used as a
 * registration between two volumes, taking into account the source
 * and destination voxel transforms.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: kteich $
 *    $Date: 2007/10/15 20:41:47 $
 *    $Revision: 1.9 $
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


#ifndef ScubaTransform_h
#define ScubaTransform_h

#include "string_fixed.h"

#include "Broadcaster.h"
#include "DebugReporter.h"
#include "IDTracker.h"
#include "Point3.h"
#include "TclCommandManager.h"
#include "Transform44.h"

class VolumeCollection;

// Static listener class so we can respond to singleton messages.
class ScubaTransformStaticTclListener : public DebugReporter,
					public TclCommandListener {

public:
  ScubaTransformStaticTclListener ();
  ~ScubaTransformStaticTclListener ();

  // Respond to Tcl commands.
  virtual TclCommandResult
    DoListenToTclCommand ( char* isCommand, int iArgc, char** iasArgv );
};

class ScubaTransform : public Transform44,
		       public IDTracker<ScubaTransform>,
		       public TclCommandListener, 
		       public Broadcaster // transformChanged
{

  friend class ScubaTransformTester;

public:

  ScubaTransform();
  virtual ~ScubaTransform();

  // Respond to Tcl commands.
  virtual TclCommandResult
    DoListenToTclCommand ( char* isCommand, int iArgc, char** iasArgv );

  // The label for this object.
  void SetLabel ( std::string const& isLabel );
  std::string const& GetLabel() const;

  // Return the main matrix values as a string.
  std::string GetValuesAsString() const;

  // Copy the transform to this one.
  ScubaTransform& operator= ( ScubaTransform const& iTransform );

  // Make this transform usable as a registration between two volumes.
  void TreatAsRegistration ( VolumeCollection const& iSourceVolume,
                             VolumeCollection const& iDestVolume );

  // Make this a nomral (non-registration) transform.
  void TreatAsNative ();

  // Return whether this is a registration transform.
  bool IsRegistration () const {
    return mIsRegistration;
  }

  // If this is a registration, return the IDs of the source and
  // destination volumes.
  int GetRegistrationSource () const;
  int GetRegistrationDestination () const;

protected:

  // Called when our values change.
  virtual void ValuesChanged ();

  // Static pointer to our static listener.
  static ScubaTransformStaticTclListener mStaticListener;

  // Our label.
  std::string msLabel;

  // For making this a registration transform. Save the source and
  // dest volume IDs so we can unmake it.
  bool mIsRegistration;
  VolumeCollection const* mSourceVolume;
  VolumeCollection const* mDestVolume;
};

std::ostream& operator << ( std::ostream&, ScubaTransform const& iTransform  );

#endif

