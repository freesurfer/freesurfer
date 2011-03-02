/**
 * @file  Path.h
 * @brief A path of Point3 objects
 *
 * This class provides a path made of Point3 objects. The path can be
 * constructed in segments, marking the end of the segment so the
 * current segment can be cleared and filled in as needed. Note: This
 * is not the best designed class (along with PathManager) so try to
 * roll with it.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:36 $
 *    $Revision: 1.7 $
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

  // Accessors.
  int GetNumVertices () const;
  Point3<T> const& GetVertexAtIndex ( int in ) const;
  Point3<T> const& GetPointAtEndOfLastSegment () const;

  // Selected flag is used by the client.
  void SetSelected ( bool ibSelected );
  bool IsSelected () const;

  // Add a vertex to the path. Notifies listeners of change.
  void AddVertex ( Point3<T> const& i );

  // Clear the entire path.
  void Clear ();

  // Segments are used by the UI for drawing paths; the last segment
  // is usually in flux, and is constantly cleared and remade as the
  // user draws the path. Then the user finalizes that segment and
  // moves onto the next one.
  void ClearLastSegment ();
  void MarkEndOfSegment ();

  // Move all points in the path by a certain amount.
  void Move ( Point3<T> const& iDelta );

  // Read and write path data to a stream.
  void ReadFromStream ( std::istream& iStream );
  void WriteToStream  ( std::ostream& ioStream ) const;

protected:

  // Broadcasts a pathChanged message when something is changed.
  void PathChanged ();

  // Broadcasts a pathVertexAdded message when a vert is added.
  void PathVertexAdded ();

  // Selected flag used by client.
  bool mbSelected;

  // Our vector of points.
  std::vector<Point3<T> > mVertices;

  // Marker for the end of the last segment.
  int mIndexOfSegmentEnd;
};

#endif
