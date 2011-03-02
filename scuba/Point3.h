/**
 * @file  Point3.h
 * @brief Simple class for a three component tuple
 *
 * This is a simple class representing a three component tuple. It has
 * the basic settors and accessors along with a stream output
 * function. All functions are defined here in the header to avoid
 * having to explicitly declare the templated instantiations.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:36 $
 *    $Revision: 1.11 $
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


#ifndef Point3_h
#define Point3_h

#include <iostream>

template <typename T>
class Point3 {

public:

  Point3 () {}

  Point3 ( T const iX, T const iY, T const iZ ) {
    m[0] = iX;
    m[1] = iY;
    m[2] = iZ;
  }

  Point3 ( T const iXYZ[3] ) {
    m[0] = iXYZ[0];
    m[1] = iXYZ[1];
    m[2] = iXYZ[2];
  }

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

  // Operators.
  inline bool operator== ( Point3<T> const& i ) {
    return ((m[0] == i.x()) && (m[1] == i.y()) && (m[2] == i.z()));
  }
  inline bool operator!= ( Point3<T> const& i ) {
    return ((m[0] != i.x()) || (m[1] != i.y()) || (m[2] != i.z()));
  }

  // Storage.
  T m[3];
};

template <typename T>
std::ostream& operator << ( std::ostream& os, Point3<T> const& iPoint ) {
  os << "(" << iPoint.x() << "," << iPoint.y() << "," << iPoint.z() << ")";
  return os;
}

#endif
