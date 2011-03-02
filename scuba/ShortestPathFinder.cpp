/**
 * @file  ShortestPathFinder.cpp
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
 *    $Revision: 1.12 $
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

#include <stdexcept>
#include <queue>
#include <vector>
#include <limits>

#include "ShortestPathFinder.h"

using namespace std;

ShortestPathFinder::ShortestPathFinder() : 
  mEdgeBias( 0.5 ),
  mzX( 0 ),
  mzY( 0 ),
  mbValidCosts( false ),
  maTotalCost( NULL ) {

  mStartPoint.Set( 0, 0 );
} 

ShortestPathFinder::~ShortestPathFinder() {

  delete maTotalCost;
}

void
ShortestPathFinder::SetDimensions ( int izX, int izY ) {

  // If these are not our current dimensions...
  if( mzX != izX || mzY != izY ) {
    
    // Save the new ones.
    mzX = izX;
    mzY = izY;
    
    // Recreate our cost array.
    delete maTotalCost;
    maTotalCost = new Array2<float>( mzX, mzY, numeric_limits<int>::max() );

    // Need to find the costs again.
    mbValidCosts = false;
  }
}

void
ShortestPathFinder::SetEdgeBias ( float iBias ) {

  if( iBias < 0 || iBias > 1.0 )
    throw runtime_error ( "Invalid edge bias; should be 0-1." );

  mEdgeBias = iBias;
}

void
ShortestPathFinder::SetStartPoint ( Point2<int> const& iStartPoint ) {

  // If this is not the same start point...
  if( mStartPoint != iStartPoint ) {

    if( iStartPoint.x() < 0 || iStartPoint.x() >= mzX ||
	iStartPoint.y() < 0 || iStartPoint.y() >= mzY )
      throw runtime_error( "Invalid start point, out of bounds." );
    
    // Set the point and invalidate the costs.
    mStartPoint = iStartPoint;
    mbValidCosts = false;
  }
}

void
ShortestPathFinder::FindPathToEndPoint ( Point2<int> const& iEndPoint,
					 std::list<Point2<int> >& ioPoints ) {

  if( iEndPoint.x() < 0 || iEndPoint.x() >= mzX ||
      iEndPoint.y() < 0 || iEndPoint.y() >= mzY )
      throw runtime_error( "Invalid start point, out of bounds." );

  // Find the costs if we haven't already.
  if( !mbValidCosts )
    FindCosts();

  // Start at the end and follow the direction table back to the
  // beginning, adding points on the way. The next step in the
  // shortest path is the neighbor with the lowest cost.
  ioPoints.clear();
  Point2<int> current( iEndPoint );
  while ( current.x() != mStartPoint.x() ||
          current.y() != mStartPoint.y() ) {

    // Push the current point to the result path.
    ioPoints.push_front( current );

    // Search for the lowest cost of its neighbors.
    float lowestCost = numeric_limits<float>::max();
    int lowestDirection = 0;
    for ( int direction = 1; direction < 9; direction++ ) {
      Point2<int> neighbor( current.x() + GetDirectionOffset(direction)[0],
			    current.y() + GetDirectionOffset(direction)[1]);
      if ( neighbor.x() >= 0 && neighbor.x() < mzX &&
	   neighbor.y() >= 0 && neighbor.y() < mzY ) {
	float cost = maTotalCost->Get( neighbor );
	if( cost < lowestCost ) {
	  lowestCost = cost;
	  lowestDirection = direction;
	}
      }
    }

    // Update current to be the neighbor with the lowest cost.
    current.Set( current.x() + GetDirectionOffset(lowestDirection)[0],
		 current.y() + GetDirectionOffset(lowestDirection)[1]);
  }
}

int const*
ShortestPathFinder::GetDirectionOffset ( int inDirection )  const {
 
  static int saDirectionOffsets[9][2] = { { 0, 0},
					  {-1, 0}, {-1, 1}, {0, 1}, { 1, 1},
					  { 1, 0}, { 1,-1}, {0,-1}, {-1,-1} };
  return saDirectionOffsets[inDirection];
}

float
ShortestPathFinder::GetDirectionFactor ( int inDirection ) const {

  static float saFactor[9] = { 0, 1, 1.4, 1, 1.4, 1, 1.4, 1, 1.4 };
  return saFactor[inDirection];
}

void
ShortestPathFinder::FindCosts () {

  // The pair is the location in the grid and the cumulative cost from
  // the start point to that location.
  priority_queue<Entry, vector<Entry>, CompareEntriesGT > q;
  
  // Initialize all costs to max.
  maTotalCost->SetAll( numeric_limits<int>::max() );

  // Set initial cost to zero and insert it in the queue.
  maTotalCost->Set( mStartPoint, 0 );
  q.push( Entry( mStartPoint, 0 ) );

  // While we have points the check...
  while ( !q.empty() ) {

    // Get the nearest entry.
    Entry currentEntry = q.top();
    q.pop();

    // Pull out the current point and cost.
    Point2<int>& current = currentEntry.first;
    float cost = currentEntry.second;

    // If this cost is lower than we total we currently have for this
    // point...
    if( cost <= maTotalCost->Get( current ) ) {

      // This is our new shortest path element. Check our neighbors to
      // update their costs.
      for ( int direction = 1; direction < 9; direction++ ) {
	Point2<int> neighbor( current.x() + GetDirectionOffset(direction)[0],
			      current.y() + GetDirectionOffset(direction)[1]);
	if ( neighbor.x() >= 0 && neighbor.x() < mzX &&
	     neighbor.y() >= 0 && neighbor.y() < mzY ) {

	  // Get the cost of this edge, multiplying the edge cost by
	  // the edge bias and the direction factor by the rest.
	  float newCost = cost +
	    ( (this->GetEdgeCost( neighbor ) * mEdgeBias) +
	      (GetDirectionFactor(direction) * (1.0 - mEdgeBias)) );

	  // If path from current point shorter than old path...
	  if ( newCost < maTotalCost->Get( neighbor ) ) {

	    // Save the new low cost.
	    maTotalCost->Set( neighbor, newCost );
	    
	    // Push it to the queue with the new cost.
	    q.push( Entry( neighbor, newCost ) );
	  }
	}
      }
    }
  }

  // We now have our costs for this start point.
  mbValidCosts = true;
}

bool
ShortestPathFinder::CompareEntriesGT::operator() ( Entry const& iA, 
						   Entry const& iB ) const {

  // Just compare the costs.
  return iA.second > iB.second;
}
