#ifndef DataCollection_h
#define DataCollection_h


#include <map>
#include <list>
#include <string>
#include "DebugReporter.h"


class DataCollection;

class DataCollection : public DebugReporter {

 public:
  typedef int ID;

  // Assigns a CollectionID and adds it to the list.
  DataCollection( std::string isLabel );

  // Removes it from the list.
  virtual ~DataCollection(); 

  // Used to poll for any displayable data at the given point.
  virtual void GetInfoAtRAS( float const iX, float const iY, float const iZ,
			     std::list<std::string> olLabels,
			     std::list<std::string> olValues ) const;

  ID GetID() const { return mID; }
  std::string GetLabel() const { return msLabel; }
  void SetLabel( std::string const isLabel ) { msLabel = isLabel; }

  // Get a DataCollection by ID.
  static DataCollection& GetDataCollection( ID const iID );

 protected:
  std::string msLabel;

  // For managing DataCollection::ID.
  ID mID;
  static ID mNextID;
  static ID GetNextID () { return mNextID++; }

  typedef std::map<DataCollection::ID,DataCollection*> CollectionIDMap;
  static CollectionIDMap mCollectionIDs;
  
};



#endif
