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
 *    $Date: 2011/03/02 00:04:11 $
 *    $Revision: 1.7 $
 *
 * Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
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


#ifndef TOPOLOGY_PATCHDISK_H
#define TOPOLOGY_PATCHDISK_H

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
