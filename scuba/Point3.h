/**
 * @file  Point3.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 02:09:14 $
 *    $Revision: 1.7 $
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
  Point3 () {}
  Point3 ( T iX, T iY, T iZ ) {
    m[0] = iX;
    m[1] = iY;
    m[2] = iZ;
  }
  Point3 ( T const iXYZ[3] ) {
    m[0] = iXYZ[0];
    m[1] = iXYZ[1];
    m[2] = iXYZ[2];
  }
  void Set ( T iX, T iY, T iZ ) {
    m[0] = iX;
    m[1] = iY;
    m[2] = iZ;
  }
  void Set ( T iXYZ[3] ) {
    m[0] = iXYZ[0];
    m[1] = iXYZ[1];
    m[2] = iXYZ[2];
  }
  void Set ( Point3<T>& i ) {
    m[0] = i.x();
    m[1] = i.y();
    m[2] = i.z();
  }
  void SetX ( T i ) {
    m[0] = i;
  }
  void SetY ( T i ) {
    m[1] = i;
  }
  void SetZ ( T i ) {
    m[2] = i;
  }
  T& operator[]( int i ) {
    return m[i];
  }
  T* xyz() {
    return m;
  }
  T x() {
    return m[0];
  }
  T y() {
    return m[1];
  }
  T z() {
    return m[2];
  }
  T m[3];
};

template <typename T>
std::ostream& operator << ( std::ostream&, Point3<T> iPoint  );


#endif
