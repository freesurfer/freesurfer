//
// DataManager.h
//
// Warning: Do not edit the following four lines.  CVS maintains them.
// Revision Author: $Author: kteich $
// Revision Date  : $Date: 2004/04/28 21:53:48 $
// Revision       : $Revision: 1.6 $

#ifndef DataManager_h
#define DataManager_h

#include <stdlib.h>
#include <list>
#include <map>
#include "string_fixed.h"
#include "DebugReporter.h"


extern "C" { 
#include "mri.h"
#include "mrisurf.h"
}


// This is the generic Loader object. It provides basic algorithm for
// keeping references to loaded data of this type, and virtual functions
// for actually loading the data and freeing it.
template <typename T> 
class DataLoader : public DebugReporter { 
 public:

  // Returns an instance to the data for the given file name, or
  // throws an error if it can't be loaded.
  T GetData( std::string const& ifnData );

  // An object should call this function when it is done with the
  // data. It may be freed, or not if it's in use by other objects.
  void ReleaseData( T* ioData );

  // Returns the number of data of this type actually loaded. (Not
  // references to a specific instance of data.
  int CountLoaded() const { return mlData.size(); }

  // Returns the number of references to this data.
  int CountReferences( T iData );

 protected:
  DataLoader();

  // Override this function to load a specfic type of data given a
  // file name. Should throw an error if it's not loadable.
  virtual T LoadData( std::string& ifnData ) = 0;

  // Override this function to free a specific type of data. Should
  // set the ioData parameter to NULL if freed.
  virtual void FreeData( T* ioData ) = 0;

  // Comparison function for determining if a specific data instance
  // matches what would be loaded for this filename.
  virtual bool DoesFileNameMatchObject( T iData, std::string& ifnData )  = 0;

  // List of data and number of references for each data.
  static std::list<T> mlData;
  static std::map<T,int> mRefs;
};


class MRILoader : public DataLoader<MRI*> {
 protected:
  MRI* LoadData( std::string& ifnData );
  void FreeData( MRI** ioMRI ) ;
  bool DoesFileNameMatchObject( MRI* iData, std::string& ifnData );
};

class MRISLoader : public DataLoader<MRIS*> {
 protected:
  MRIS* LoadData( std::string& ifnData );
  void FreeData( MRIS** ioMRI ) ;
  bool DoesFileNameMatchObject( MRIS* iData, std::string& ifnData );
};



class DataManager : public DebugReporter {

 public: 

  static DataManager& GetManager();

  static MRILoader& GetMRILoader();
  static MRISLoader& GetMRISLoader();

 protected:
  DataManager();
  
  //  static MRILoader mMRILoader;
  //  static MRISLoader mMRISLoader;
};




#endif

