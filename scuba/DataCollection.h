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
  DataLocation () { 
    mRAS[0] = 0;
    mRAS[1] = 0;
    mRAS[2] = 0;
  }
  DataLocation ( float const iRAS[3] ) {
    mRAS[0] = iRAS[0]; 
    mRAS[1] = iRAS[1]; 
    mRAS[2] = iRAS[2]; 
  }
  DataLocation ( const DataLocation& iLoc ) {
    mRAS[0] = iLoc.RAS(0);
    mRAS[1] = iLoc.RAS(1);
    mRAS[2] = iLoc.RAS(2);
  }
  ~DataLocation () {}
  float* RAS() { return mRAS; }
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

  DataCollection();
  virtual ~DataCollection(); 

  // If the normal DataLocation is not enough, should subclass to
  // create specific DataLocation. Basically for caching RAS -> data
  // index, be it MRI coords or vertex index or whatever.
  virtual DataLocation& MakeLocationFromRAS ( float const iRAS[3] );

  // used to poll for any displayable data at the given point.
  virtual void GetInfo( DataLocation& iLoc,
			std::map<std::string,std::string>& iLabelValues );

  // Should return a type description unique to the subclass.
  virtual std::string GetTypeDescription() { return "BaseCollection"; }

  std::string GetLabel() const { return msLabel; }
  void SetLabel( std::string const isLabel );
  
  // Return the bounds of the data in RAS coords. 0=xmin, 1=xmax,
  // 2=ymin, etc.
  virtual void GetDataRASBounds ( float oBounds[6] );

  virtual TclCommandResult
    DoListenToTclCommand ( char* isCommand, int iArgc, char** iasArgv );

  // Handle broadcast messages.
  virtual void
    DoListenToMessage ( std::string isMessage, void* iData );

  // Get a list of ROI IDs that belong to this data collection.
  std::vector<int> GetROIList ();
  int GetNumberOfROIs () { return mROIMap.size(); }
  bool IsROIInThisCollection ( int iROIID );
  
  // Create a new ROI and assign it to this collection. Return its ID.
  int NewROI ();

  // Tell this collection to select this ROI.
  void SelectROI ( int iROIID );
  int GetSelectedROI () { return mSelectedROIID; }
  void DeleteROI ( int iROIID );

  // Called by NewROI, should be subclassed to return specific ROI type.
  virtual ScubaROI* DoNewROI ();

  // The Data <-> World transform.
  virtual void SetDataToWorldTransform ( int iTransformID );
  int GetDataToWorldTransform ();

  // Returns a best guess value increment for a GUI.
  virtual float GetPreferredValueIncrement ();

  // Suppresses the dataChanged message. Use when changing a lot of
  // voxels in a row that don't need updates in between. Will call
  // DataChanged() at the end.
  void BeginBatchChanges ();
  void EndBatchChanges ();

  // Passes batch change messages down to the current ROI.
  void BeginBatchROIChanges ();
  void EndBatchROIChanges ();

protected:
  std::string msLabel;

  int mSelectedROIID;
  std::map<int,ScubaROI*> mROIMap;
  
  // For self to call when data has changed.
  virtual void DataChanged ();
  bool mbSuspendDataChangedMessage;
 
  // The data to world transform. Should be applied to all requests
  // for data at RAS points.
  ScubaTransform* mDataToWorldTransform;

};



#endif
