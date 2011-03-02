/**
 * @file  Volume3.cpp
 * @brief Simple class for a 3D volume of elements
 *
 * Template class for creating 3D volumes of a data type with settors
 * and accessors.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:39 $
 *    $Revision: 1.10 $
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


#include "string_fixed.h"
#include <stdexcept>

#ifdef __cplusplus
extern "C"
{
#endif
#include <stdlib.h> // calloc
#ifdef __cplusplus
}
#endif

#include "Volume3.h"

using namespace std;

template <typename T>
Volume3<T>::Volume3 ( int izX, int izY, int izZ, T iInitialValue ) {

  mzX = izX;
  mzY = izY;
  mzZ = izZ;

  // Allocate our storage.
  mData = (T***) calloc( mzZ, sizeof(T**) );
  for ( int nZ = 0; nZ < mzZ; nZ++ ) {
    mData[nZ] = (T**) calloc( mzY, sizeof(T*) );
    for ( int nY = 0; nY < mzY; nY++ ) {
      mData[nZ][nY] = (T*) calloc( mzX, sizeof(T) );
      for ( int nX = 0; nX < mzX; nX++ ) {
        mData[nZ][nY][nX] = iInitialValue;
      }
    }
  }

}



template <typename T>
Volume3<T>::~Volume3 () {

  if ( NULL != mData ) {
    for ( int nZ = 0; nZ < mzZ; nZ++ ) {
      for ( int nY = 0; nY < mzY; nY++ ) {
        free( mData[nZ][nY] );
      }
      free( mData[nZ] );
    }
    free( mData );
  }
  mData = NULL;
}

template <typename T>
T
Volume3<T>::Get ( int inX, int inY, int inZ ) const {

  if ( inX >= 0 && inX < mzX &&
       inY >= 0 && inY < mzY &&
       inY >= 0 && inZ < mzZ ) {
    return mData[inZ][inY][inX];
  } else {
    throw runtime_error( "Out of bounds" );
  }
}

template <typename T>
T
Volume3<T>::Get_Unsafe ( int inX, int inY, int inZ ) const {

  return mData[inZ][inY][inX];
}

template <typename T>
void
Volume3<T>::GetBounds( int& ozX, int& ozY, int& ozZ ) const {

  ozX = mzX;
  ozY = mzY;
  ozZ = mzZ;
}


template <typename T>
void
Volume3<T>::SetAll ( T const iValue ) {

  for ( int nZ = 0; nZ < mzZ; nZ++ ) {
    for ( int nY = 0; nY < mzY; nY++ ) {
      for ( int nX = 0; nX < mzX; nX++ ) {
        mData[nZ][nY][nX] = iValue;
      }
    }
  }
}

template <typename T>
void
Volume3<T>::Set ( int inX, int inY, int inZ, T const iValue ) {

  if ( inX >= 0 && inX < mzX &&
       inY >= 0 && inY < mzY &&
       inY >= 0 && inZ < mzZ ) {
    mData[inZ][inY][inX] = iValue;
  } else {
    throw runtime_error( "Out of bounds" );
  }
}

template <typename T>
void
Volume3<T>::Set_Unsafe ( int inX, int inY, int inZ, T const iValue ) {

  mData[inZ][inY][inX] = iValue;
}


