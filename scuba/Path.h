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
    return mVertices[mIndexOfSegmentEnd]; }

  T GetSquaredDistanceOfClosestPoint ( T iMinDistance, Point3<T>& iWhere );

  void Move ( Point3<T>& iDelta );

  int GetNumVertices () { return mVertices.size(); }
  Point3<T>& GetVertexAtIndex ( int in ) { return mVertices[in]; }

  // Broadcasts a pathChanged message when something is changed.
  void PathChanged ();

  // Broadcasts a pathVertexAdded message when a vert is added.
  void PathVertexAdded ();

 protected:
  std::vector<Point3<T> > mVertices;
  int mIndexOfSegmentEnd;
};

#endif
