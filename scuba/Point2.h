/**
 * @file  Point2.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 02:09:14 $
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


#ifndef Point2_h
#define Point2_h

#include <iostream>


template <typename T>
class Point2 {
public:
  Point2 () {}
  Point2 ( T iX, T iY ) {
    m[0] = iX;
    m[1] = iY;
  }
  Point2 ( T const iXY[2] ) {
    m[0] = iXY[0];
    m[1] = iXY[1];
  }
  void Set ( T iX, T iY ) {
    m[0] = iX;
    m[1] = iY;
  }
  void Set ( T iXY[2] ) {
    m[0] = iXY[0];
    m[1] = iXY[1];
  }
  void Set ( Point2<T>& i ) {
    m[0] = i.x();
    m[1] = i.y();
  }
  T& operator[]( int i ) {
    return m[i];
  }
  T* xy() {
    return m;
  }
  T x () {
    return m[0];
  }
  T y () {
    return m[1];
  }
  T m[2];
};


template <typename T>
std::ostream& operator << ( std::ostream&, Point2<T> iPoint  );

#endif
