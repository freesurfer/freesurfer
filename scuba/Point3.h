/**
 * @file  Point3.h
 * @brief Simple class for a three component tuple
 *
 * This is a simple class representing a three component tuple. It has
 * the basic settors and accessors along with a stream output function.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: kteich $
 *    $Date: 2007/10/10 18:43:36 $
 *    $Revision: 1.8 $
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


#ifndef Point3_h
#define Point3_h

#include <iostream>

template <typename T>
class Point3 {

public:

  Point3 ();
  Point3 ( T const iX, T const iY, T const iZ );
  Point3 ( T const iXYZ[3] );

  // Settors.
  inline void Set ( T const iX, T const iY, T const iZ ) {
    m[0] = iX;
    m[1] = iY;
    m[2] = iZ;
  }
  inline void Set ( T const iXYZ[3] ) {
    m[0] = iXYZ[0];
    m[1] = iXYZ[1];
    m[2] = iXYZ[2];
  }
  inline void Set ( Point3<T> const& i ) {
    m[0] = i.x();
    m[1] = i.y();
    m[2] = i.z();
  }
  inline void SetX ( T const i ) {
    m[0] = i;
  }
  inline void SetY ( T const i ) {
    m[1] = i;
  }
  inline void SetZ ( T const i ) {
    m[2] = i;
  }
  inline T& operator[]( int i ) {
    return m[i];
  }
  inline T* xyz() {
    return m;
  }

  // Accessors.
  inline const T& operator[]( int i ) const {
    return m[i];
  }
  inline T const* xyz() const {
    return m;
  }
  inline T x() const {
    return m[0];
  }
  inline T y() const {
    return m[1];
  }
  inline T z() const {
    return m[2];
  }

  // Storage.
  T m[3];
};

template <typename T>
std::ostream& operator << ( std::ostream&, Point3<T> const& iPoint  );


#endif
