/**
 * @file  ShortestPathFinder.cpp
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 02:09:15 $
 *    $Revision: 1.9 $
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


#include "ShortestPathFinder.h"

//  Portions (c) Copyright 2005 Brigham and Women's Hospital (BWH) All
//  Rights Reserved.
//  See docs/license.slicer
//  or http://www.slicer.org/copyright/copyright.txt for details.

using namespace std;

ShortestPathFinder::ShortestPathFinder() {

  mStraightBias = 1;
  mEdgeBias = 0.5;
  mLongestEdge = 0;
  mzX = mzY = 0;
  mQueue = NULL;
  maCost = NULL;
  maDir = NULL;
  maDone = NULL;
  mDebug = false;
}

ShortestPathFinder::~ShortestPathFinder() {

  if ( NULL != mQueue ) {
    delete mQueue;
  }
  if ( NULL != maCost ) {
    delete maCost;
  }
  if ( NULL != maDir ) {
    delete maDir;
  }
  if ( NULL != maDone ) {
    delete maDone;
  }
}

void
ShortestPathFinder::SetDimensions ( int izX, int izY, int iLongestEdge ) {

  mzX = izX;
  mzY = izY;
  mLongestEdge = iLongestEdge;

  if ( NULL != mQueue ) {
    delete mQueue;
  }
  if ( NULL != maCost ) {
    delete maCost;
  }
  if ( NULL != maDir )  {
    delete maDir;
  }
  if ( NULL != maDone ) {
    delete maDone;
  }

  mQueue = new circularQueue( mzX, mzY, iLongestEdge );
  maCost = new Array2<float>( mzX, mzY, iLongestEdge * mzX * mzY ); // max cost
  maDir  = new Array2<int>( mzX, mzY, 0 );
  maDone = new Array2<bool>( mzX, mzY, false );
}


void
ShortestPathFinder::FindPath ( Point2<int>& iStartPoint,
                               Point2<int>& iEndPoint,
                               std::list<Point2<int> >& ioPoints ) {

  int aDirectionOffsets[9][2] = { { 0, 0},
                                  {-1, 0}, {-1, 1}, {0, 1}, { 1, 1},
                                  { 1, 0}, { 1,-1}, {0,-1}, {-1,-1} };
  float aFactor[9] = { 0, 1, 1.4, 1, 1.4, 1, 1.4, 1, 1.4 };

  // Set begin and ending points.
  Point2<int> begin( iStartPoint );
  Point2<int> end( iEndPoint );

  // Set initial cost to zero and insert it in the queue.
  maCost->Set( begin, 0 );
  mQueue->Insert( begin, 0 );

  // While end point not done, keep checking out neighbors of current point
  float currentCost = 0;
  while ( !maDone->Get( end ) ) {

    // Get min vertex from Q
    listElement* min = mQueue->GetListElement( (int)currentCost );
    if ( NULL == min ) {
      cerr << "GetListElement returned NULL" << endl;
      return;
    }
    Point2<int> current( min->Coord );

    // Update the current path cost
    currentCost = maCost->Get( current );

    // Mark this point done.
    maDone->Set( current, true );

    // remove it from Q
    mQueue->Remove( min );

    // Check out neighbors.
    for ( int direction = 1; direction < 9; direction++ ) {

      Point2<int> neighbor( current.x() + aDirectionOffsets[direction][0],
                            current.y() + aDirectionOffsets[direction][1]);

      if ( neighbor.x() >= 0 && neighbor.x() < mzX &&
           neighbor.y() >= 0 && neighbor.y() < mzY ) {

        // Get the cost of this edge = value of weight volume *
        // factor.  We add 0.001 to the edge cost because if it's
        // 0, there's no preference for straight lines, since
        // diagonal lines will have the same cost, which in some
        // cases makes a really weird looking line.
        float newCost = currentCost +
                        ( (this->GetEdgeCost( neighbor ) * mEdgeBias + 0.001) *
                          aFactor[direction] * (1.0 - mEdgeBias) );

        // If path from current point shorter than old path
        if ( newCost < maCost->Get( neighbor ) ) {

          // lower the cumulative cost to reach this neighbor
          maCost->Set( neighbor, newCost );

          // store new short path direction
          maDir->Set( neighbor, direction );

          // remove this neighbor from Q if in it
          mQueue->Remove( neighbor );

          // then put it in the proper place in Q for its new cost
          mQueue->Insert( neighbor, (int)newCost );
        }
      }
    }

#if 0
    if ( mDebug ) {
      cerr << "after pulling" << current << endl;
      for ( int nY = 0; nY < mzY; nY++ ) {
        for ( int nX = 0; nX < mzX; nX++ ) {
          if ( nX == current.x() && nY == current.y() ) {
            cerr.width( 4 );
            cerr.fill( '-' );
            cerr << maCost->Get( nX, nY ) << "-";
          } else {
            cerr.width( 4 );
            cerr.fill( ' ' );
            cerr << maCost->Get( nX, nY ) << " ";
          }
        }
        cerr << endl;
      }
      cerr << endl;
    }
#endif
  }

  // Start at the end and follow the direction table back to the
  // beginning, adding points on the way.
  ioPoints.clear();
  Point2<int> current( end );
  Point2<int> next;
  while ( current.x() != begin.x() ||
          current.y() != begin.y() ) {

    ioPoints.push_front( current );

    next.Set( current.x() - aDirectionOffsets[maDir->Get(current)][0],
              current.y() - aDirectionOffsets[maDir->Get(current)][1]);

    current.Set( next.x(), next.y() );
  }

#if 0
  if ( mDebug ) {
    cerr << "final cost array and path" << endl;
    for ( int nY = 0; nY < mzY; nY++ ) {
      for ( int nX = 0; nX < mzX; nX++ ) {

        bool bInPath = false;
        list<Point2<int> >::iterator tPoint;
        for ( tPoint = ioPoints.begin(); tPoint != ioPoints.end();
              ++tPoint ) {
          if ( nX == (*tPoint).x() && nY == (*tPoint).y() ) {
            bInPath = true;
            break;
          }
        }
        if ( bInPath ) {
          cerr.width( 7 );
          cerr.precision( 4 );
          cerr.fill( '-' );
          cerr << maCost->Get( nX, nY ) << "-";
        } else {
          cerr.width( 7 );
          cerr.precision( 4 );
          cerr.fill( ' ' );
          cerr << maCost->Get( nX, nY ) << " ";
        }
      }
      cerr << endl;
    }
    cerr << endl;
  }
#endif
}


linkedList::linkedList(int x, int y)
    : Array2<listElement>(x,y) {
  for (int i = 0; i < x; i++) {
    for (int j = 0; j < y; j++) {
      this->Element(i,j)->Coord.Set( i, j );
    }
  }
}


circularQueue::circularQueue(int x, int y, int buckets) {
  this->A = new linkedList(x,y);
  this->C = buckets;
  this->Circle = new listElement[this->C+1];
  // link each bucket into its circle
  for (int i=0; i<C+1; i++) {
    this->Circle[i].Prev = this->Circle[i].Next = &this->Circle[i];
  }
};

circularQueue::~circularQueue() {
  if (this->A) delete this->A;
  if (this->Circle) delete[] this->Circle;
};

void
circularQueue::Insert ( Point2<int>& iLocation, int iCost ) {

  int nBucket = this->GetBucket( iCost );

  listElement *el = this->A->Element( iLocation );
  // insert el at the top of the list from the bucket
  el->Next = this->Circle[nBucket].Next;
  if (el->Next == NULL) {
    cout << "ERROR. bucket is NULL, not linked to self." << endl;
  }
  this->Circle[nBucket].Next->Prev = el;
  this->Circle[nBucket].Next = el;
  el->Prev = &this->Circle[nBucket];
}

void
circularQueue::Remove ( Point2<int>& iLocation ) {

  listElement *el = this->A->Element( iLocation );
  this->Remove( el );
}

void
circularQueue::Remove( listElement *el ) {

  // if el is in linked list
  if (el->Prev != NULL) {
    if (el->Next == NULL) {
      cout <<"ERROR. el->Next is NULL."<< endl;
      return;
    }

    el->Next->Prev = el->Prev;
    el->Prev->Next = el->Next;

    // clear el's pointers
    el->Prev = el->Next = NULL;
  }

  return;
}

listElement*
circularQueue::GetListElement ( int iCost ) {

  int nBucket = FindMinBucket( iCost );

  // return the last one in the linked list.
  if (this->Circle[nBucket].Prev == NULL) {
    cout << "ERROR. Unlinked list." << endl;
    return NULL;
  }
  if (this->Circle[nBucket].Next == &this->Circle[nBucket]) {
    cout << "ERROR. Empty linked list." << endl;
    return NULL;
  }

  return this->Circle[nBucket].Prev;
}

int
circularQueue::GetBucket( int iCost ) {

  if ( iCost < 0 ) {
    cout << "ERROR: negative cost of " << iCost << endl;
  }

  // return remainder
  return div(iCost,this->C+1).rem;
}

int
circularQueue::FindMinBucket( int iCost ) {

  int nBucket = this->GetBucket( iCost );
  int count = 0;

  int cost = iCost;
  while ( this->Circle[nBucket].Next == &this->Circle[nBucket] &&
          count <= this->C ) {

    // search around the Q for the next vertex
    cost++;
    nBucket = this->GetBucket( cost );
    count++;
  }

  // have we looped all the way around?
  if (count > this->C) {
    cout << "ERROR. Empty Q." << endl;
  }

  if (this->Circle[nBucket].Prev == &this->Circle[nBucket]) {
    cout <<"ERROR. Prev not linked to bucket." << endl;
  }

  return nBucket;
}
