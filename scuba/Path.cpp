/**
 * @file  Path.cpp
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
 *    $Author: kteich $
 *    $Date: 2007/10/16 20:19:22 $
 *    $Revision: 1.11 $
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

#include <limits>

#include "Path.h"
#include "VectorOps.h"

using namespace std;

template <typename T>
Path<T>::Path () :
  Broadcaster( "Path" ),
  mbSelected( false ),
  mIndexOfSegmentEnd( 0 ) {
}

template <typename T>
int
Path<T>::GetNumVertices () const {

  return mVertices.size();
}

template <typename T>
Point3<T> const& 
Path<T>::GetVertexAtIndex ( int in ) const {
  
  return mVertices[in];
}

template <typename T>
Point3<T> const& 
Path<T>::GetPointAtEndOfLastSegment () const {
  
  return mVertices[mIndexOfSegmentEnd];
}

template <typename T>
void
Path<T>::SetSelected ( bool ibSelected ) {
  
  mbSelected = ibSelected;
}

template <typename T>
bool
Path<T>::IsSelected () const {

  return mbSelected;
}

template <typename T>
void
Path<T>::AddVertex ( Point3<T> const& i ) {
  
  // Add this point to our vector.
  mVertices.push_back( i );

  // Let us notify listeners.
  PathVertexAdded();
}

template <typename T>
void
Path<T>::Clear () {

  // Clear our vector.
  mVertices.clear();

  // Let us notify listeners.
  PathChanged();
}


template <typename T>
void
Path<T>::ClearLastSegment () {

  // Pop entries from the index of the last segment to the end.
  typename vector<Point3<T> >::iterator tEndOfLastSegment;
  tEndOfLastSegment = mVertices.begin();
  tEndOfLastSegment += mIndexOfSegmentEnd+1;
  mVertices.erase( tEndOfLastSegment, mVertices.end() );

  // Let us notify listeners.
  PathChanged();
}

template <typename T>
void
Path<T>::MarkEndOfSegment () {

  // Take the current size-1.
  mIndexOfSegmentEnd = mVertices.size() - 1;
}

template <typename T>
void
Path<T>::Move ( Point3<T> const& iDelta ) {

  typename vector<Point3<T> >::iterator tPoint;
  for ( tPoint = mVertices.begin(); tPoint != mVertices.end(); ++tPoint ) {
    Point3<float>& point = (*tPoint);
    point.Set( point.x() + iDelta.x(),
               point.y() + iDelta.y(),
               point.z() + iDelta.z() );
  }

  // Notify listeners.
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

  // Read the number of vertices.
  int cVertices;
  iStream >> cVertices;

  // Read values and construct points.
  for ( int nVertex = 0; nVertex < cVertices; nVertex++ ) {
    T x, y, z;
    iStream >> x >> y >> z;
    Point3<T> point( x, y, z );
    AddVertex( point );
  }

  // Mark this as the end of the segment.
  MarkEndOfSegment();
}

template <typename T>
void
Path<T>::WriteToStream ( ostream& ioStream ) const {

  // Write the number of vertices.
  ioStream << mVertices.size() << endl;

  // Write all the values to the stream.
  typename vector<Point3<T> >::const_iterator tPoint;
  for ( tPoint = mVertices.begin(); tPoint != mVertices.end(); ++tPoint ) {
    Point3<float> const& point = (*tPoint);
    ioStream << point.x() << " " << point.y() << " " << point.z() << endl;
  }
}
