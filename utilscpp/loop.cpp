/**
 * @file  loop.cpp
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 02:09:19 $
 *    $Revision: 1.2 $
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


#include "topology/loop.h"

Loop::Loop(void) {
  npoints=0;
  maxpoints=0;
  points=0;
}

Loop::Loop(int maxpts) {
  points=0;
  Alloc(maxpts);
}

Loop::~Loop(void) {
  if (points) delete [] points;
}

void Loop::Alloc(int maxpts) {
  if (points) delete [] points;
  points = 0 ;
  maxpoints = maxpts;
  points = new int[maxpoints];
  npoints=0;
}

#define INCREASE_FACTOR 1.2
void Loop::_ReAlloc(int maxpts) {
  if (maxpts < maxpoints)
    maxpoints = int(maxpoints*INCREASE_FACTOR+1);
  else maxpoints = maxpts;

  int *new_points = new int[maxpoints];
  memcpy(new_points,points,npoints*sizeof(int));
  if (points) delete [] points;
  points=new_points;
}

void Loop::AddPoint(int pt) {
  if (npoints==maxpoints) {
    if (maxpoints==0) _ReAlloc(10);
    else _ReAlloc();
  }
  points[npoints++]=pt;
}


void Loop::Print() const {
  for (int n = 0 ; n < npoints-1 ; n++)
    cout << points[n] << "->";
  cout << points[npoints-1]<< "." << endl;
}

int Loop::End() {
  return points[npoints-1];
}

void Loop::Pop() {
  if (npoints==0) return;
  npoints--;
}

int Loop::Replace(int pt, int new_pt) {
  int nreplaced=0;
  for (int n = 0 ; n < npoints ; n++)
    if (points[n]==pt) {
      points[n]=new_pt;
      nreplaced++;
    }
  return nreplaced;
}

const Loop & Loop::operator=(const Loop& loop) {
  if (maxpoints < loop.npoints) _ReAlloc(loop.npoints);
  npoints = loop.npoints;
  for (int n = 0 ; n < npoints ; n++)
    points[n] = loop.points[n];
  return loop;
}
