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

  // Returns whether a point is in this path. Path must be an inplane
  // polygon and be closed, otherwise the result won't be useful.
  bool PointInPath ( Point3<T>& iPoint );
  
  void Move ( Point3<T>& iDelta );

  int GetNumVertices () { return mVertices.size(); }
  Point3<T>& GetVertexAtIndex ( int in ) { return mVertices[in]; }

  // Broadcasts a pathChanged message when something is changed.
  void PathChanged ();

  // Broadcasts a pathVertexAdded message when a vert is added.
  void PathVertexAdded ();

  void ReadFromStream ( std::istream& iStream );
  void WriteToStream  ( std::ostream& ioStream );

  void SetSelected ( bool ibSelected ) { mbSelected = ibSelected; }
  bool IsSelected () { return mbSelected; }

 protected:
  bool mbSelected;

  std::vector<Point3<T> > mVertices;
  int mIndexOfSegmentEnd;
};

#endif
