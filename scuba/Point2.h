/**
 * @file  Point2.h
 * @brief Simple class for a two component tuple
 *
 * This is a simple class representing a two component tuple. It has
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


#ifndef Point2_h
#define Point2_h

#include <iostream>

template <typename T>
class Point2 {

public:
  
  Point2 () {}

  Point2 ( T const iX, T const iY ) {
    m[0] = iX;
    m[1] = iY;
  }

  Point2 ( T const iXY[2] ) {

    m[0] = iXY[0];
    m[1] = iXY[1];
  }

  // Settors.
  inline void Set ( T const iX, T const iY ) {
    m[0] = iX;
    m[1] = iY;
  }
  inline void Set ( T const iXY[2] ) {
    m[0] = iXY[0];
    m[1] = iXY[1];
  }
  inline void Set ( Point2<T> const& i ) {
    m[0] = i.x();
    m[1] = i.y();
  }
  inline T& operator[]( int i ) {
    return m[i];
  }
  inline T* xy() {
    return m;
  }

  // Accessors.
  inline T const& operator[]( int i ) const {
    return m[i];
  }
  inline T const* xy() const {
    return m;
  }
  inline T x () const {
    return m[0];
  }
  inline T y () const {
    return m[1];
  }

  // Operators.
  inline bool operator== ( Point2<T> const& i ) {
    return ((m[0] == i.x()) && (m[1] == i.y()));
  }
  inline bool operator!= ( Point2<T> const& i ) {
    return ((m[0] != i.x()) || (m[1] != i.y()));
  }

  // Storage.
  T m[2];
};


template <typename T>
std::ostream& operator << ( std::ostream& os, Point2<T> const& iPoint  ) {
  os << "(" << iPoint.x() << "," << iPoint.y() << ")";
  return os;
}

#endif
