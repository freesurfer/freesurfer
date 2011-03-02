/**
 * @file  Utilities.h
 * @brief Collection of utilities
 *
 * Some utilities that don't fit in anywhere else.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:39 $
 *    $Revision: 1.8 $
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


#ifndef Utilities_h
#define Utilities_h

#include <list>
#include <vector>
#include "Point2.h"
#include "Point3.h"


class Utilities {

public:

  // Return a list of points in the 2D coord plane that are on the
  // line from A to B.
  static void FindPointsOnLine2d ( int const iPointA[2],
				   int const iPointB[2],
                                   std::list<Point2<int> >& oPoints );

  // Find the distance from a point to the closest point on a line or
  // segment.
  static float DistanceFromLineToPoint3f ( Point3<float> const& iLineA,
					   Point3<float> const& iLineB,
					   Point3<float> const& iPoint );
  static float DistanceFromSegmentToPoint3f ( Point3<float> const& iLineA,
					      Point3<float> const& iLineB,
					      Point3<float> const& iPoint );

  // Splits a string into individual strings based on a delimiter, and
  // returns the count of delimiters found (so the oResults will have
  // that +1 elements in it.)
  // Following is written by Paul J. Weiss,
  // http://www.codeproject.com/string/stringsplit.asp
  static int SplitString( std::string const& isInput,
                          std::string const& isDelimiter,
                          std::vector<std::string>& oResults );

};


#endif
