/**
 * @file  ShortestPathFinder.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 02:09:15 $
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


#ifndef ShortestPathFinder_h
#define ShortestPathFinder_h

//  Portions (c) Copyright 2005 Brigham and Women's Hospital (BWH) All
//  Rights Reserved.
//  See docs/license.slicer
//  or http://www.slicer.org/copyright/copyright.txt for details.

#include <list>
#include "Point2.h"
#include "Array2.h"
#include "DebugReporter.h"

// This class calculates shortest paths in a 2d plane given a cost to
// go to each point. (This is possibly incorrectly called 'edge cost',
// but is actually the cost to get to that point from any other.) The
// path from point A to B is returned as A+1 to B, not including the
// original point A.
// To use it, subclass Shortest and redefine GetEdgeCost if you want.


class listElement {
public:
  listElement *Prev;
  listElement *Next;
  Point2<int> Coord;
  listElement() {
    this->Prev=NULL;
    this->Next=NULL;
  };
};


class linkedList : public Array2<listElement> {
public:
  linkedList ( int izX, int izY );
};

class circularQueue {
public:
  circularQueue ( int izX, int izY, int icBuckets );
  ~circularQueue();

  void Insert ( Point2<int>& iLocation, int iCost );
  void Remove ( Point2<int>& iLocation );
  void Remove ( listElement *el );
  listElement *GetListElement ( int iCost );

private:
  int GetBucket ( int iCost );
  int FindMinBucket ( int iCost );

  linkedList *A;
  listElement *Circle;
  int C;
};



class ShortestPathFinder : public DebugReporter {

public:

  ShortestPathFinder();
  virtual ~ShortestPathFinder();

  void SetDimensions ( int izX, int izY, int iLongestEdge );

  void SetStraightBias ( float iBias ) {
    mStraightBias = iBias;
  }
  void SetEdgeBias ( float iBias ) {
    mEdgeBias = iBias;
  }

  void FindPath ( Point2<int>& iStartPoint, Point2<int>& iEndPoint,
                  std::list<Point2<int> >& ioPoints );

  virtual float GetEdgeCost ( Point2<int>& ) {
    return 1.0;
  }

  void SetDebug ( bool iDebug ) {
    mDebug = iDebug;
  }


protected:
  float mStraightBias;
  float mEdgeBias;
  int mzX, mzY;
  float mLongestEdge;
  circularQueue *mQueue;
  Array2<float> *maCost;
  Array2<int> *maDir;
  Array2<bool> *maDone;
  bool mDebug;
};


#endif
