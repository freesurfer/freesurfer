//
// DataManager.h
//
// Warning: Do not edit the following four lines.  CVS maintains them.
// Revision Author: $Author: kteich $
// Revision Date  : $Date: 2003/10/06 00:52:03 $
// Revision       : $Revision: 1.1 $

#ifndef DataManager_h
#define DataManager_h

#include <stdlib.h>
#include <list>
#include <map>
#include <string>

extern "C" { 
#include "mri.h"
#include "mrisurf.h"
}


template <class T> 
class Loader { 
 public:

  T GetData( char const* iFileName );
  void ReleaseData( T* ioData );
  int CountLoaded() const { return mlData.size(); }

 protected:
  virtual T LoadData( std::string& ifnData ) {};
  virtual void FreeData( T* ioData ) {};
  virtual bool DoesFileNameMatchObject( T iData, std::string& ifnData ) {};

  std::list<T> mlData;
  std::map<T,int> mRefs;
};



class MRILoader : public Loader<MRI*> {
 protected:

  MRI* LoadData( std::string& ifnData ) { 
    char* fnMRI = strdup( ifnData.c_str() );
    MRI* mri = MRIread( fnMRI ); 
    free( fnMRI );
    if( NULL == mri ) {
      throw (char const*) "Couldn't load MRI.";
    }
    return mri;
  }

  void FreeData( MRI** ioMRI ) { 
    MRIfree( ioMRI ); 
  }

  virtual bool DoesFileNameMatchObject( MRI* iData, std::string& ifnData ) {
    std::string fnCur( iData->fname );
    return (fnCur == ifnData);
  }
};


class DataManager {

 public: 

  static DataManager& GetManager();

  MRI* GetMRI( char const* ifnMRI );
  void ReleaseMRI( MRI** ioMRI );
  int CountLoadedMRIs() const;

 protected:
  DataManager();
  
  std::list<MRI*> mlMRI;
  std::map<MRI*,int> mMRIRefs ;

  Loader<MRI*> mMRILoader;

};




#endif
