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
  DataLocation ( float const iRAS[3] ) {
    mRAS[0] = iRAS[0]; 
    mRAS[1] = iRAS[1]; 
    mRAS[2] = iRAS[2]; 
  }
  ~DataLocation () {}
  float* RAS() { return mRAS; }
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

  virtual DataLocation& MakeLocationFromRAS ( float const iRAS[3] );

  // used to poll for any displayable data at the given point.
  virtual void GetInfo( DataLocation& iLoc,
			std::map<std::string,std::string>& iLabelValues );

  // Should return a type description unique to the subclass.
  virtual std::string GetTypeDescription() { return "BaseCollection"; }

  std::string GetLabel() const { return msLabel; }
  void SetLabel( std::string const isLabel ) { msLabel = isLabel; }
  
  virtual TclCommandResult
    DoListenToTclCommand ( char* isCommand, int iArgc, char** iasArgv );

  // Handle broadcast messages.
  virtual void
    DoListenToMessage ( std::string isMessage, void* iData );

  std::vector<int> GetROIList ();

  int NewROI ();
  void SelectROI ( int iROIID );
  virtual ScubaROI* DoNewROI ();

  int GetSelectedROI () { return mSelectedROIID; }

  virtual void SetDataToWorldTransform ( int iTransformID );
  int GetDataToWorldTransform ();

  // Suppresses the dataChanged message. Use when changing a lot of
  // voxels in a row that don't need updates in between. Will call
  // DataChanged() at the end.
  void BeginBatchChanges ();
  void EndBatchChanges ();

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
