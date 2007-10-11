/**
 * @file  Array2.h
 * @brief Simple class for a 2D array
 *
 * This is a simple class representing a 2D array, accessible by x,y
 * pairs or by a Point2<int>.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: kteich $
 *    $Date: 2007/10/11 21:45:43 $
 *    $Revision: 1.3 $
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

  Array2 ( int izX, int izY );
  Array2 ( int izX, int izY, T const iInitValue );
  ~Array2 ();

  // Settors
  void SetAll ( T const iValue );
  inline void Set ( int inX, int inY, T const iValue ) {
    this->mArray[ArrayIndexToInternalIndex(inX,inY)] = iValue;
  }
  inline void Set ( Point2<int> const& iLocation, T const iValue ) {
    this->mArray[ArrayIndexToInternalIndex(iLocation)] = iValue;
  }
  inline T& operator () ( int const inX, int const inY ) {
    return this->mArray[ArrayIndexToInternalIndex(inX,inY)];
  }
  inline T& operator () ( Point2<int> const& iLocation ) {
    return this->mArray[ArrayIndexToInternalIndex(iLocation)];
  }

  // Accessors
  inline T const& Get ( int const inX, int const inY ) const {
    return this->mArray[ArrayIndexToInternalIndex(inX,inY)];
  }
  inline T& Get ( int const inX, int const inY ) {
    return this->mArray[ArrayIndexToInternalIndex(inX,inY)];
  }
  inline T const& Get ( Point2<int> const& iLocation ) const {
    return this->mArray[ArrayIndexToInternalIndex(iLocation)];
  }
  inline T& Get ( Point2<int> const& iLocation ) {
    return this->mArray[ArrayIndexToInternalIndex(iLocation)];
  }
  inline T const& operator () ( int const inX, int const inY ) const {
    return this->mArray[ArrayIndexToInternalIndex(inX,inY)];
  }
  inline T const& operator () ( Point2<int> const& iLocation ) const {
    return this->mArray[ArrayIndexToInternalIndex(iLocation)];
  }

private:

  // Convert x,y or point coords to an internal array index.
  inline int ArrayIndexToInternalIndex ( int iX, int iY ) const {
    return iX + (iY * this->mzX);
  }
  inline int ArrayIndexToInternalIndex ( Point2<int> const& i ) const {
    return i.x() + (i.y() * this->mzX);
  }

  // Storage
  T* mArray;
  int mzX, mzY;
};


#endif
