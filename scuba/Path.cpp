#include "Path.h"
#include "VectorOps.h"

using namespace std;

template <typename T>
Path<T>::Path () {

  mIndexOfSegmentEnd = 0;
  mbSelected = false;
}

template <typename T>
void 
Path<T>::Clear () {

  mVertices.clear();

  PathChanged();
}

template <typename T>
void 
Path<T>::AddVertex ( Point3<T>& i ) {

  mVertices.push_back( i );

  PathVertexAdded();
}

template <typename T>
void
Path<T>::MarkEndOfSegment () {

  mIndexOfSegmentEnd = mVertices.size() - 1;
}

template <typename T>
void
Path<T>::ClearLastSegment () {

  for( int cToPop = mVertices.size() - mIndexOfSegmentEnd-1;
       cToPop > 0; cToPop-- ) {
    
    mVertices.pop_back();
  }

  PathChanged();
}

template <typename T>
T
Path<T>::GetSquaredDistanceOfClosestPoint ( T iMinDistance, 
					    Point3<T>& iWhere ) {
  
  T minDistance = iMinDistance;

  typename vector<Point3<T> >::iterator tPoint;
  for( tPoint = mVertices.begin(); tPoint != mVertices.end(); ++tPoint ) {
    Point3<float>& point = (*tPoint);
    T distance = 
      (iWhere.x() - point.x()) * (iWhere.x() - point.x()) +
      (iWhere.y() - point.y()) * (iWhere.y() - point.y()) +
      (iWhere.z() - point.z()) * (iWhere.z() - point.z());
    if( distance < minDistance ) {
      minDistance = distance;
    }
  }
  
  return minDistance;
}

template <typename T>
bool
Path<T>::PointInPath ( Point3<T>& iPoint ) {

  // Make sure we have at least 3 points.
  if( mVertices.size() < 3 ) 
    throw runtime_error( "Path too short." );

  // Make sure we're all on the same plane and find that plane.
  bool bSamePlane[3] = {true, true, true};
  Point3<T> min;
  min.Set( 10000, 10000, 10000 );
  typename vector<Point3<T> >::iterator tPoint;
  for( tPoint = mVertices.begin(); tPoint != mVertices.end(); ++tPoint ) {
    Point3<T>& first = mVertices[0];
    Point3<T>& point = (*tPoint);
    if( point[0] != first[0] ) bSamePlane[0] = false;
    if( point[1] != first[1] ) bSamePlane[1] = false;
    if( point[2] != first[2] ) bSamePlane[2] = false;

    // Also find a min.
    if( point[0] < min[0] ) min[0] = point[0];
    if( point[1] < min[1] ) min[1] = point[1];
    if( point[2] < min[2] ) min[2] = point[2];
  }

  // If we have no same planes, return false.
  if( !bSamePlane[0] && !bSamePlane[1] && !bSamePlane[2] ) {
    throw runtime_error( "Not all in same plane." );
  }

  // Make sure we're closed.
  Point3<T>& first = mVertices[0];
  Point3<T>& last = mVertices[mVertices.size()-1];
  if( first != last ) 
    throw runtime_error( "Not closed." );

  // Make sure input point is in same plane.
  if( bSamePlane[0] && first[0] != iPoint[0] ||
      bSamePlane[1] && first[1] != iPoint[1] ||
      bSamePlane[2] && first[2] != iPoint[2] ) {
    throw runtime_error( "Input point not in right plane." );
  }
  
  // Find a point on the outside edge.
  Point3<T> p2 = min;

  // Now count intersections.
  Point3<T> x;
  int cIntersections = 0;
  int cVertices = mVertices.size();
  for( int nVertex = 0; nVertex < cVertices-1; nVertex++ ) {
    Point3<T> q1 = mVertices[nVertex];
    Point3<T> q2 = mVertices[nVertex+1];
    
    VectorOps::IntersectionResult rIntersect = 
      VectorOps::SegmentIntersectsSegment( iPoint, p2, q1, q2, x );
    if( VectorOps::intersect == rIntersect ) {
      cIntersections++;
    }
  }

  if( cIntersections % 2 == 0 ) 
    return false;
  else
    return true;
}

template <typename T>
void 
Path<T>::Move ( Point3<T>& iDelta ) {

  typename vector<Point3<T> >::iterator tPoint;
  for( tPoint = mVertices.begin(); tPoint != mVertices.end(); ++tPoint ) {
    Point3<float>& point = (*tPoint);
    point.Set( point.x() + iDelta.x(), 
	       point.y() + iDelta.y(), 
	       point.z() + iDelta.z() );
  }
  
  PathChanged();
}

template <typename T>
void 
Path<T>::PathChanged () {
  
  // Broadcast this change.
  int id = this->GetID();
  SendBroadcast( "pathChanged", (void*)&id );
}

template <typename T>
void 
Path<T>::PathVertexAdded () {
  
  // Broadcast this change.
  int id = this->GetID();
  SendBroadcast( "pathVertexAdded", (void*)&id );
}

template <typename T>
void
Path<T>::ReadFromStream ( istream& iStream ) {

  int cVertices;
  iStream >> cVertices;

  for( int nVertex = 0; nVertex < cVertices; nVertex++ ) {
    T x, y, z;
    iStream >> x >> y >> z;
    Point3<T> point( x, y, z );
    AddVertex( point );
  }

  MarkEndOfSegment();
}

template <typename T>
void
Path<T>::WriteToStream ( ostream& ioStream ) {

  ioStream << mVertices.size() << endl;

  typename vector<Point3<T> >::iterator tPoint;
  for( tPoint = mVertices.begin(); tPoint != mVertices.end(); ++tPoint ) {
    Point3<float>& point = (*tPoint);
    ioStream << point.x() << " " << point.y() << " " << point.z() << endl;
  }
}
