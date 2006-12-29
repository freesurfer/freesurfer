/**
 * @file  segment.cpp
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 02:09:19 $
 *    $Revision: 1.3 $
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


#include "topology/segment.h"

Segment::Segment(void) {
  npoints=0;
  marked=0;
  maxpoints=NUMBER_OF_POINTS;
  points = new int[maxpoints];
}

Segment::~Segment(void) {
  if (points) delete [] points;
}

int Segment::GetMark() const {
  return marked;
}
void Segment::SetMark(int m) {
  marked=m;
}

void Segment::_ReallocSegment(int new_maxpoints) {
  if (new_maxpoints < 0 || new_maxpoints < maxpoints)
    new_maxpoints = int(maxpoints*INCREASE_NUMBER_OF_POINTS+1);

  int *new_points = new int[new_maxpoints];
  for (int i = 0 ; i < maxpoints ; i++)
    new_points[i]=points[i];

  maxpoints=new_maxpoints;
  delete [] points;
  points=new_points;
}

void Segment::AddPoint(int pt) {
  if (npoints == maxpoints)
    _ReallocSegment();
  points[npoints++]=pt;
}

int Segment::size() const {
  return npoints;
}
void Segment::clear() {
  npoints=0;
}

void Segment::AddSegment(Segment *s) {
  int new_maxpoints = npoints+s->npoints;
  if (new_maxpoints > maxpoints) _ReallocSegment(new_maxpoints);

  //copy the segment s
  memcpy(&points[npoints],s->points,s->npoints*sizeof(int));
  npoints += s->npoints;

  s->clear();
}

const int* Segment::GetPointList(void) const {
  return points;
}

void Segment::Transfer(Segment &b) {
  npoints=b.npoints;
  maxpoints=b.maxpoints;
  delete [] points;
  points=b.points;
  b.maxpoints=0;
  b.npoints=0;
  b.points=NULL;
}
