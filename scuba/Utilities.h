/**
 * @file  Utilities.h
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


#ifndef Utilities_h
#define Utilities_h

#include <list>
#include <vector>
#include "Point2.h"
#include "Point3.h"


class Utilities {

public:

  static void FindPointsOnLine2d ( int iPointA[2], int iPointB[2],
                                   int iThickness,
                                   std::list<Point2<int> >& oPoints );

  static void FindPointsOnLine3d ( int iPointA[3], int iPointB[3],
                                   std::list<Point3<int> >& oPoints );
  static void FindPointsOnLine3f ( float iPointA[3], float iPointB[3],
                                   std::list<Point3<float> >& oPoints );

  static float DistanceFromLineToPoint3f ( Point3<float>& iLineA,
      Point3<float>& iLineB,
      Point3<float>& iPoint );
  static float DistanceFromSegmentToPoint3f ( Point3<float>& iLineA,
      Point3<float>& iLineB,
      Point3<float>& iPoint );

  // Following is written by
  // Paul J. Weiss, http://www.codeproject.com/string/stringsplit.asp
  static int SplitString( const std::string& isInput,
                          const std::string& isDelimiter,
                          std::vector<std::string>& oResults );

};


#endif
