/**
 * @file  DataCollection.h
 * @brief Base abstract data collection object
 *
 * This class is the base class for all data objects in Scuba. A 'collection'
 * is usually a primary data object, such as a volume or surface, and its
 * associated data, such as ROIs. This file also defines the DataLocation
 * class, which is encapsulates different coordinate spaces such as RAS and
 * data indicies.

 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:36 $
 *    $Revision: 1.26 $
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


#ifndef DataCollection_h
#define DataCollection_h


#include <map>
#include <vector>
#include "string_fixed.h"
#include "DebugReporter.h"
#include "IDTracker.h"
#include "TclCommandManager.h"
#include "ScubaROI.h"
#include "ScubaTransform.h"


// This class is used by clients to reference a sample point. It is a
// way of avoiding access functions for multiple data types (i.e. one
// for RAS, one for index, etc) and allowing the DataCollection to
// cache coordinate conversions.
class DataLocation {
  friend class DataCollectionTester;
  friend class DataCollection;
public:

  // Ctors.
  DataLocation ();
  DataLocation ( float const iRAS[3] );
  DataLocation ( DataLocation const& iLoc );
  ~DataLocation () {}

  // Accessors.
  float const* RAS() const { return mRAS; }
  float RAS ( int in ) const { return mRAS[in]; }

protected:
  float mRAS[3];
};

class DataCollection : public DebugReporter,
      public IDTracker<DataCollection>,
      public TclCommandListener,
      public Listener,    // transformChanged
      public Broadcaster  // dataChanged
{

  friend class DataCollectionTester;

public:

  // The Ctor sets the mDataToWorldTransform to a default identity
  // transform with ID 0; if not found, it creates it.
  DataCollection();
  virtual ~DataCollection();

  // If the normal DataLocation is not enough, should subclass to
  // create specific DataLocation. Basically for caching RAS -> data
  // index, be it MRI coords or vertex index or whatever.
  virtual DataLocation MakeLocationFromRAS ( float const iRAS[3] ) const;

  // Should return a type description unique to the subclass.
  virtual std::string GetTypeDescription() const {
    return "BaseCollection";
  }

  // Get and set our human-readable label.
  std::string const& GetLabel() const;
  void SetLabel( std::string const& isLabel );

  // Return the bounds of the data in RAS coords. 0=xmin, 1=xmax,
  // 2=ymin, etc.
  virtual void GetDataRASBounds ( float oBounds[6] ) const;

  // Respond to Tcl commands.
  virtual TclCommandResult
    DoListenToTclCommand ( char* isCommand, int iArgc, char** iasArgv );

  // Handle broadcast messages.
  virtual void
    DoListenToMessage ( std::string isMessage, void* iData );
  
  // ROI management. Each DataCollection has a set of ROIs that belong
  // to it, and each knows how to make the right kind of ROI for it.

  // Create a new ROI and assign it to this collection. Return its
  // ID. This calls the internal DoNewROI() function that subclasses
  // can override to make their own ROI types. We take over ownership
  // of the ROI pointer.
  int NewROI ();

  // Delete a given ROI.
  void DeleteROI ( int iROIID );

  // Get a list of ROI IDs that belong to this data collection.
  std::vector<int> GetROIList () const;
  int GetNumberOfROIs () const;
  bool IsROIInThisCollection ( int iROIID ) const;

  // A collection can have a selected ROI. This lets it perform
  // certain actions on the selected ROI.
  void SelectROI ( int iROIID );
  int GetSelectedROI () const;

  // A DataCollection subclass will most likely define a custom space
  // for whatever coordinates are in its data fields, such as
  // column/row/slice space for a volume. "Data" space is an
  // intermediate space that allows a user-transform (i.e. a display
  // transform) to be inserted between this space and the final world
  // space. "World" space is equivalent to RAS space, since in Scuba,
  // everything that is drawn on the screen is meant to be in RAS
  // space. This superclass only deals with the Data <-> World
  // transform.

  // Set the ID of the transform (this should be the ID of a
  // ScubaTransform).
  virtual void SetDataToWorldTransform ( int iTransformID );

  // Return the ID of the transform
  int GetDataToWorldTransform () const;


  // Returns a best guess value increment for a GUI. This is used by
  // editor tools that change this data type's values. This is
  // convenience function and may not be applicable to all subclasses.
  virtual float GetPreferredValueIncrement () const;


  // Suppresses the dataChanged message. Use when changing a lot of
  // voxels in a row that don't need updates in between. Will call
  // DataChanged() at the end.
  void BeginBatchChanges ();
  void EndBatchChanges ();

  // Passes batch change messages down to the current ROI.
  void BeginBatchROIChanges ();
  void EndBatchROIChanges ();

protected:

  // Called by NewROI, should be subclassed to return specific ROI type.
  virtual ScubaROI* DoNewROI ();

  // Our human-readable label.
  std::string msLabel;

  // Our map of ROIs that belong to us. We own the pointers in this map.
  std::map<int,ScubaROI*> mROIMap;

  // Our selected ROI.
  int mSelectedROIID;

  // For self to call when data has changed.
  virtual void DataChanged ();

  // Don't call DataChanged if true (set by the
  // {Begin|End}BatchChanges calls);
  bool mbSuspendDataChangedMessage;

  // The data to world transform. Should be applied to all requests
  // for data at RAS points.
  ScubaTransform* mDataToWorldTransform;

};



#endif
