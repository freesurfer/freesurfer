#include <iostream>
#include <string>
#include <list>
#include "DataManager.h"

using namespace std;


template <typename T>
T DataLoader<T>::GetData( string const& ifnData ) {

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

template <typename T>
void
DataLoader<T>::ReleaseData( T* ioData ) {

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


MRILoader DataManager::mMRILoader;
MRISLoader DataManager::mMRISLoader;

DataManager::DataManager() {
}

DataManager& 
DataManager::GetManager() {
  static DataManager sManager;
  return sManager;
}

// This line necessary to generate the right code.
template DataLoader<MRI*>;

MRI* 
MRILoader::LoadData( std::string& ifnData ) { 

  // Use MRIread to load the MRI object. Need to make a non-const, 
  // c-string copy of the file name.
  char* fnMRI = strdup( ifnData.c_str() );
  MRI* mri = MRIread( fnMRI ); 
  free( fnMRI );
  if( NULL == mri ) {
    throw (char const*) "Couldn't load MRI.";
  }
  return mri;
}

void
MRILoader::FreeData( MRI** ioMRI ) { 

  // Call MRIfree. This will set *ioMRI to NULL if successful.
  MRIfree( ioMRI ); 
}

bool
MRILoader::DoesFileNameMatchObject( MRI* iData, std::string& ifnData ) {

  // Look at the fname member in the MRI structure and compare it
  // to the file name we're getting.
  std::string fnCur( iData->fname );
  return (fnCur == ifnData);
}


// This line necessary to generate the right code.
template DataLoader<MRIS*>;

MRIS* 
MRISLoader::LoadData( std::string& ifnData ) { 

  // Use MRISread to load the MRIS object. Need to make a non-const, 
  // c-string copy of the file name.
  char* fnMRIS = strdup( ifnData.c_str() );
  MRIS* mris = MRISread( fnMRIS ); 
  free( fnMRIS );
  if( NULL == mris ) {
    throw (char const*) "Couldn't load MRIS.";
  }
  return mris;
}

void
MRISLoader::FreeData( MRIS** ioMRIS ) { 

  // Call MRISfree. This will set *ioMRIS to NULL if successful.
  MRISfree( ioMRIS ); 
}

bool
MRISLoader::DoesFileNameMatchObject( MRIS* iData, std::string& ifnData ) {

  // Look at the fname member in the MRI structure and compare it
  // to the file name we're getting.
  std::string fnCur( iData->fname );
  return (fnCur == ifnData);
}
