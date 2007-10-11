/**
 * @file  Point2.h
 * @brief Simple class for a two component tuple
 *
 * This is a simple class representing a two component tuple. It has
 * the basic settors and accessors along with a stream output function.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: kteich $
 *    $Date: 2007/10/11 21:45:44 $
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


#ifndef Point2_h
#define Point2_h

#include <iostream>

template <typename T>
class Point2 {

public:
  
  Point2 ();
  Point2 ( T const iX, T const iY );
  Point2 ( T const iXY[2] );

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
std::ostream& operator << ( std::ostream&, Point2<T> const& iPoint  );

#endif
