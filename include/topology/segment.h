/**
 * @file  segment.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 02:09:02 $
 *    $Revision: 1.4 $
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


#ifndef TOPOLOGY_SEGMENT_H
#define TOPOLOGY_SEGMENT_H

#ifdef __cplusplus

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

#endif
