/**
 * @file  patchdisk.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 02:09:02 $
 *    $Revision: 1.6 $
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


#ifndef TOPOLOGY_PATCHDISK_H
#define TOPOLOGY_PATCHDISK_H

#ifdef __cplusplus

#define MAX_EXTRA_VERTICES 100 //66
#define MAX_EXTRA_FACES 150 //128


#include "surface.h"

class PatchDisk
{
  void _Init();
  void _Alloc(int which_patch);
public:
  Surface disk;
  Loop ring,init_ring;
  int *vtrans ;
  int *ftrans ;

  PatchDisk();
  PatchDisk(int which_patch);
  PatchDisk(const string s):disk(s),init_ring(10)
  {
    vtrans = new int[disk.nvertices];
    ftrans = new int[disk.nfaces];
    _Init();
  };
  ~PatchDisk(void);

  void Init();
  void Create(int which_patch);
};


#endif

#endif
