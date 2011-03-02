/**
 * @file  ShortestPathFinder.h
 * @brief Calculates shortest path in 2D
 *
 * This class calculates shortest paths in a 2D plane given a cost to
 * go to each point. An edge bias factor allows difference importance
 * to be placed on the edge cost vs. the distance (1.0 for cardinal
 * points, 1.4 for diagonals). The path from point A to B is returned
 * as A+1 to B, not including the original point A.  To use it,
 * subclass Shortest and redefine GetEdgeCost if you want.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:39 $
 *    $Revision: 1.10 $
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

//  Portions (c) Copyright 2005 Brigham and Women's Hospital (BWH) All
//  Rights Reserved.
//  See docs/license.slicer
//  or http://www.slicer.org/copyright/copyright.txt for details.

#ifndef ShortestPathFinder_h
#define ShortestPathFinder_h

#include <list>

#include "Array2.h"
#include "DebugReporter.h"
#include "Point2.h"

class ShortestPathFinder : public DebugReporter {

public:

  ShortestPathFinder();
  virtual ~ShortestPathFinder();

  // Set the dimensions of the range. If these are already the current
  // dimensions, does nothing.
  void SetDimensions ( int izX, int izY );

  // Set the edge bias, from 0 to 1. This determines how much weight
  // to give to the edge cost vs the direction cost (which is 1 for
  // cardinal transitions and 1.4 for diagonal transitions.
  void SetEdgeBias ( float iBias );

  // Set the start point and find all the paths. This performs the
  // searching algorithm. After this is called, you can call
  // FindPathToEndPoint to use this data. If this function is called
  // to change the start point, the paths must be recalculated.
  void SetStartPoint ( Point2<int> const& iStartPoint );

  // Find the shortest path between the previously set start and end
  // point. ioPoints will be cleared and the result path will be
  // stored in it as a list of points from the point AFTER the start
  // point up to and INCLUDING the end point.
  void FindPathToEndPoint ( Point2<int> const& iEndPoint,
			    std::list<Point2<int> >& ioPoints );

  // Subclasses should override this to get the edge cost for a
  // vertex.
  virtual float GetEdgeCost ( Point2<int> const& ) const {
    return 1.0;
  }

protected:

  // Find the costs from the start point.
  void FindCosts ();

  // Returns a direction offset array of x,y. or a distance factor for
  // each direction: for cardinal directions and 1.4 for
  // diagonals. Valid numbers are from 1 to 8, representing 8
  // directions you can travel. No input checking.
  int const* GetDirectionOffset ( int inDirection ) const;
  float GetDirectionFactor ( int inDirection ) const;

  // Type used by the priority_queue in FindPath.
  typedef std::pair<Point2<int>,float> Entry;

  // Comparator used by the priority_queue in FindPath.
  class CompareEntriesGT {
  public: 
    bool operator() ( Entry const&, Entry const& ) const;
  };

  // Edge bias, from 0 to 1.
  float mEdgeBias;

  // Dimensions.
  int mzX, mzY;

  // The current start point.
  Point2<int> mStartPoint;

  // Whether or not the cost array is valid.
  bool mbValidCosts;

  // We calc this in SetStartPoint to store the total shortest cost
  // from the start point to any point in the array.
  Array2<float>* maTotalCost;
};


#endif
