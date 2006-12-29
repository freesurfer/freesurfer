/**
 * @file  Volume3.cpp
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 02:09:15 $
 *    $Revision: 1.6 $
 *
 * Copyright (C) 2002-2007,
 * The General Hospital Corporation (Boston, MA). 
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 * Bug reports: analysis-bugs@nmr.mgh.harvard.edu
 *
 */


#include "string_fixed.h"
#include <stdexcept>

#include "Volume3.h"

using namespace std;

template <typename T>
Volume3<T>::Volume3 ( int izX, int izY, int izZ, T iInitialValue ) {

  mzX = izX;
  mzY = izY;
  mzZ = izZ;

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
void
Volume3<T>::SetAll ( T iValue ) {

  for ( int nZ = 0; nZ < mzZ; nZ++ ) {
    for ( int nY = 0; nY < mzY; nY++ ) {
      for ( int nX = 0; nX < mzX; nX++ ) {
        mData[nZ][nY][nX] = iValue;
      }
    }
  }
}

template <typename T>
T
Volume3<T>::Get ( int inX, int inY, int inZ ) {

  if ( inX >= 0 && inX < mzX &&
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

  if ( inX >= 0 && inX < mzX &&
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


