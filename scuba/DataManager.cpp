#include <iostream>
#include <string>
#include <list>
#include "DataManager.h"

using namespace std;




T
Loader::GetData( char const* ifnData ) {

  list<T>::iterator tData;
  string fnData( ifnData );

  for( tData = mlData.begin(); tData != mlData.end(); ++tData ) {
    T data = *tData;
    string fnCurData( ifnData );
    if( this->DoesFileNameMatchObject( data, fnCurData ) ) {
      mRefs[data]++;
      return data;
    }
  }

  T data = this->LoadData( fnData );

  mRefs[data] = 1;

  mlData.push_back( data );

  return data;
}


void
Loader::ReleaseData( T* ioData ) {

  list<T>::iterator tData;

  for( tData = mlData.begin(); tData != mlData.end(); ++tData ) {
    T data = *tData;

    if( data == *ioData ) {

      if( 1 == mRefs[data] ) {
	mlData.remove( data );
	this->FreeData( &data );
      }

      mRefs[data]--;
      *ioData = NULL;
      return;
    }
  }

  throw (char const*) "Couldn't find data";
}



DataManager::DataManager() {

}

DataManager& 
DataManager::GetManager() {

  static DataManager sManager;

  return sManager;
}

MRI*
DataManager::GetMRI( char const* ifnMRI ) {

  list<MRI*>::iterator tMRI;
  string fnMRI( ifnMRI );

#ifdef DEBUG
  cerr << "GetMRI( " << ifnMRI << " )" << endl;
#endif

  for( tMRI = mlMRI.begin(); tMRI != mlMRI.end(); ++tMRI ) {
    MRI* mri = *tMRI;
    string fnCurMRI( mri->fname );
    if( fnMRI == fnCurMRI ) {
      
#ifdef DEBUG
      cerr << "\tFound existing with " << mMRIRefs[mri] << " references,"
	   << " returning" << endl;
#endif

      mMRIRefs[mri]++;
      return mri;
    }
  }

#ifdef DEBUG
      cerr << "\tNot found, calling MRIread" << endl;
#endif

  char* fnMRI2 = strdup( ifnMRI );
  MRI* mri = MRIread( fnMRI2 );
  if( NULL == mri ) {
    throw (char const*) "Couldn't load MRI.";
  }

  mMRIRefs.insert(map<MRI*,int>::value_type(mri, 1));
  mMRIRefs[mri] = 1;

  mlMRI.push_back( mri );

  return mri;
}

void
DataManager::ReleaseMRI( MRI** ioMRI ) {

  list<MRI*>::iterator tMRI;

#ifdef DEBUG
  cerr << "ReleaseMRI()" << endl;
#endif

  for( tMRI = mlMRI.begin(); tMRI != mlMRI.end(); ++tMRI ) {
    MRI* mri = *tMRI;

    if( mri == *ioMRI ) {

#ifdef DEBUG
      cerr << "\tFound. Number of refs is " << mMRIRefs[mri] << endl;
#endif

      if( 1 == mMRIRefs[mri] ) {
	mlMRI.remove( mri );
	MRIfree( &mri );
      }

      mMRIRefs[mri]--;
      *ioMRI = NULL;
      return;
    }
  }

  throw (char const*) "Couldn't find MRI";
}

int
DataManager::CountLoadedMRIs() const {

#ifdef DEBUG
  cerr << "CountLoadedMRIs()" << endl;
  cerr << "\tmlMRI.size() = " << mlMRI.size() << endl;
#endif

  return mlMRI.size();
}



