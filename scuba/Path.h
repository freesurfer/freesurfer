/**
 * @file  Path.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 02:09:14 $
 *    $Revision: 1.5 $
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


#ifndef Path_h
#define Path_h

#include <vector>
#include "Point3.h"
#include "IDTracker.h"
#include "Broadcaster.h"

template <typename T>
class Path : public IDTracker<Path<T> >,
      public Broadcaster { // pathChanged <id>, pathVertexAdded <id>

  friend class PathTester;

public:
  Path ();

  void Clear ();

  void ClearLastSegment ();

  void MarkEndOfSegment ();

  // Notifies listeners of change.
  void AddVertex ( Point3<T>& i );

  Point3<T>& GetPointAtEndOfLastSegment () {
    return mVertices[mIndexOfSegmentEnd];
  }

  T GetSquaredDistanceOfClosestPoint ( T iMinDistance, Point3<T>& iWhere );

  // Returns whether a point is in this path. Path must be an inplane
  // polygon and be closed, otherwise the result won't be useful.
  bool PointInPath ( Point3<T>& iPoint );

  void Move ( Point3<T>& iDelta );

  int GetNumVertices () {
    return mVertices.size();
  }
  Point3<T>& GetVertexAtIndex ( int in ) {
    return mVertices[in];
  }

  // Broadcasts a pathChanged message when something is changed.
  void PathChanged ();

  // Broadcasts a pathVertexAdded message when a vert is added.
  void PathVertexAdded ();

  void ReadFromStream ( std::istream& iStream );
  void WriteToStream  ( std::ostream& ioStream );

  void SetSelected ( bool ibSelected ) {
    mbSelected = ibSelected;
  }
  bool IsSelected () {
    return mbSelected;
  }

protected:
  bool mbSelected;

  std::vector<Point3<T> > mVertices;
  int mIndexOfSegmentEnd;
};

#endif
