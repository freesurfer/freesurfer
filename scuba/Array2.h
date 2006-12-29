/**
 * @file  Array2.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 02:09:13 $
 *    $Revision: 1.2 $
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


#ifndef Array2_h
#define Array2_h

#include "Point2.h"

template <typename T>
class Array2 {
public:
  Array2 ( int const izX, int const izY );
  Array2 ( int const izX, int const izY, T const iInitValue );
  ~Array2 () {
    if ( this->mArray ) delete[] this->mArray;
  };

  // Get/Set functions
  T Get ( int const inX, int const inY ) const {
    return this->mArray[inX + inY * this->mzX];
  };
  void Set ( int const inX, int const inY, T const iValue ) {
    this->mArray[inX + inY * this->mzX] = iValue;
  };

  void SetAll ( T const iValue );

  // return the array element itself
  T& operator () ( int const inX, int const inY ) {
    return this->mArray[inX + inY * this->mzX];
  };

  // return a pointer to the array element
  T* Element ( int const inX, int const inY ) {
    return( this->mArray + inX + inY * this->mzX );
  };


  T Get ( Point2<int>& iLocation ) const {
    return this->mArray[iLocation.x() + iLocation.y() * this->mzX];
  };
  void Set ( Point2<int>& iLocation, T const iValue ) {
    this->mArray[iLocation.x() + iLocation.y() * this->mzX] = iValue;
  };

  T& operator () ( Point2<int>& iLocation ) {
    return this->mArray[ iLocation.x() + iLocation.y() * this->mzX];
  };
  T* Element ( Point2<int>& iLocation ) {
    return( this->mArray + iLocation.x() + iLocation.y() * this->mzX );
  };


private:
  T* mArray;
  int mzX, mzY;
};


#endif
