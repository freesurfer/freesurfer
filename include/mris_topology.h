/**
 * @file  mris_topology.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 02:09:00 $
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


// include guard
#ifndef MRIS_TOPOLOGY_INCLUDED_H
#define MRIS_TOPOLOGY_INCLUDED_H

// The following is usable from C
#ifdef __cplusplus
extern "C"
{
#endif

#include "mrisurf.h"

  //functions
  bool MRIScorrectDefect(MRIS *mris, int defect_number,TOPOFIX_PARMS &parms);
  bool MRIScorrectPatchTopology(MRIS* &mris,TOPOFIX_PARMS &parms);
  bool MRISincreaseEuler(MRIS* &mris,TOPOFIX_PARMS &parms);
  int MRISgetEuler(MRIS *mris, int defect_number = -1);
  int MRISgetEulerNumber(const MRIS *mris, const int *list_of_faces, int nfs);
  MRIP* MRIPextractFromMRIS(MRIS *mris, int defect_number);
  void MRISinitSurface(MRIS *mris);
  bool MRISaddMRIP(MRIS *mris_dst, MRIP *mrip);
  MRIS *MRISduplicateOver(MRIS *mris,int mode = 0);
  void MRIScopyHeader(MRIS *mris_src,MRIS *mris_dst);

#ifdef __cplusplus
}
#endif



// C++ portion starts here
#ifdef __cplusplus

extern "C"
{
#include "mrisurf.h"
#include "error.h"
}
#include "topology/surface.h"


MRIP *MRIPalloc(int nvertices, int nfaces);
void MRIPfree(MRIP **mrip);
MRIP *MRIPclone(MRIP *src);

Surface *MRIStoSurface(MRIS *mris);
MRIS * SurfaceToMRIS(Surface *surface, MRIS *mris);


#endif

#endif
