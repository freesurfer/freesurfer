#include "Path.h"

using namespace std;

template <typename T>
Path<T>::Path () {

  mIndexOfSegmentEnd = 0;
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
  SendBroadcast( "pathChanged", (void*)&mID );
}

template <typename T>
void 
Path<T>::PathVertexAdded () {
  
  // Broadcast this change.
  SendBroadcast( "pathVertexAdded", (void*)&mID );
}

