/*
 *
 * Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
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


#ifndef TOPOLOGY_SEGMENT_H
#define TOPOLOGY_SEGMENT_H

#include "globals.h"

#define NUMBER_OF_POINTS 10
#define INCREASE_NUMBER_OF_POINTS 1.2

class Segment
{
private:
  int npoints,maxpoints;
  int* points;

  int nvertices,nfaces,nedges,euler;
  int marked;

  void _ReallocSegment(int new_maxpoints=-1);

public:
  // constructor / destructor
  Segment(void);
  ~Segment(void);

  int size() const;
  void clear();
  void AddPoint(int pt);
  void AddSegment(Segment *s);
  void Transfer(Segment &b);
  int GetEuler();
  const int* GetPointList(void) const;
  int GetMark() const;
  void SetMark(int m);
};

#endif
