/**
 * @file  loop.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 02:09:02 $
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


#ifndef TOPOLOGY_LOOP_H
#define TOPOLOGY_LOOP_H

#ifdef __cplusplus

#include "globals.h"

class Loop
{
private:
  void _ReAlloc(int maxpts=-1);
public:
  int npoints,maxpoints;
  int *points;
public:
  //constructor/destructor
  Loop(void);
  Loop(int maxpts);
  ~Loop(void);

  void Alloc(int maxpts);
  void AddPoint(int pt);
  void Print() const;
  int End();
  void Pop();
  int Replace(int pt, int new_pt);
  int operator[](int n)
  {
    ASSERT((n >= 0 ) && (n < npoints));
    return points[n];
  }
  const Loop & operator=(const Loop& loop);
};

#endif

#endif
