#include <string>
#include <stdexcept>

#include "Volume3.h"

using namespace std;

template <typename T>
Volume3<T>::Volume3 ( int izX, int izY, int izZ, T iInitialValue ) {

  mzX = izX; 
  mzY = izY; 
  mzZ = izZ; 

  mData = (T***) calloc( mzZ, sizeof(T**) );
  for( int nZ = 0; nZ < mzZ; nZ++ ) {
    mData[nZ] = (T**) calloc( mzY, sizeof(T*) );
    for( int nY = 0; nY < mzY; nY++ ) {
      mData[nZ][nY] = (T*) calloc( mzX, sizeof(T) );
      for( int nX = 0; nX < mzX; nX++ ) {
	mData[nZ][nY][nX] = iInitialValue;
      }
    }
  }

}

template <typename T>
Volume3<T>::~Volume3 () {

  if( NULL != mData ) {
    for( int nZ = 0; nZ < mzZ; nZ++ ) {
      for( int nY = 0; nY < mzY; nY++ ) {
	free( mData[nZ][nY] );
      }
      free( mData[nZ] );
    }
    free( mData );
  }
}

template <typename T>
T
Volume3<T>::Get ( int inX, int inY, int inZ ) {

  if( inX >= 0 && inX < mzX &&
      inY >= 0 && inY < mzY &&
      inY >= 0 && inZ < mzZ ) {
    return mData[inZ][inY][inX];
  } else {
    throw runtime_error( "Out of bounds" );
  }
}

template <typename T>
void 
Volume3<T>::Set ( int inX, int inY, int inZ, T iValue ) {

  if( inX >= 0 && inX < mzX &&
      inY >= 0 && inY < mzY &&
      inY >= 0 && inZ < mzZ ) {
    mData[inZ][inY][inX] = iValue;
  } else {
    throw runtime_error( "Out of bounds" );
  }
}

template <typename T>
T
Volume3<T>::Get_Unsafe ( int inX, int inY, int inZ ) {

  return mData[inZ][inY][inX];
}

template <typename T>
void 
Volume3<T>::Set_Unsafe ( int inX, int inY, int inZ, T iValue ) {

  mData[inZ][inY][inX] = iValue;
}


