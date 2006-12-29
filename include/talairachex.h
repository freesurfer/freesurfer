/**
 * @file  talairachex.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 02:09:00 $
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


//
// talairachex.h
//
// Warning: Do not edit the following four lines.  CVS maintains them.
// Revision Author: $Author: nicks $
// Revision Date  : $Date: 2006/12/29 02:09:00 $
// Revision       : $Revision: 1.3 $
/////////////////////////////////////////////////////////////////////////////
// new talairach related routines with Ex
//  takes lta as the talairach transform (use LTAreadEx routine)
//  don't rely on COR volume type
/////////////////////////////////////////////////////////////////////////////
#ifndef c_talairach_h
#define c_talairach_h

#include "matrix.h"
#include "mri.h"
#include "transform.h"

MATRIX *MtalairachFromVoxel(MRI *mri_src, const LTA *lta);
MATRIX *MtalVoxelFromVoxel(MRI *mri_src, const LTA *lta);
MATRIX *MvoxelFromTalairach(MRI *mri_dst, const LTA *lta);
MATRIX *MvoxelFromTalVoxel(MRI *mri_dst, const LTA *lta);
MATRIX *MRASFromTalVoxel(MRI *mri, const LTA *lta);
void TransformWithMatrix(MATRIX *mat, Real x, Real y, Real z, Real *px, Real *py, Real *pz);

int
ModifyTalairachCRAS(MRI *mri_tal, const LTA *lta);
int
MRIvoxelToTalairachEx(MRI *mri_src, Real xv, Real yv, Real zv,
                      Real *pxt, Real *pyt, Real *pzt, const LTA *lta);
int
MRItalairachToVoxelEx(MRI *mri_dst, Real xt, Real yt, Real zt,
                      Real *pxv, Real *pyv, Real *pzv, const LTA *lta);
int
MRItalairachVoxelToWorldEx(MRI *mri_dst, Real xtv, Real ytv, Real ztv,
                           Real *pxw, Real *pyw, Real *pzw, const LTA *lta);
int
MRIvoxelToTalairachVoxelEx(MRI *mri_src, Real xv, Real yv, Real zv,
                           Real *pxt, Real *pyt, Real *pzt, const LTA *lta);
int
MRItalairachVoxelToVoxelEx(MRI *mri_dst, Real xv, Real yv, Real zv,
                           Real *pxnv, Real *pynv, Real *pznv, const LTA *lta) ;
MRI *
MRItoTalairachEx(MRI *mri_src, MRI *mri_tal, const LTA *lta);
MRI *
MRItoTalairachExInterp(MRI *mri_src, MRI *mri_tal, const LTA *lta, int interp);
MRI *
MRIfromTalairachEx(MRI *mri_tal, MRI *mri_dst, const LTA *lta);
int
MRIeraseTalairachPlaneNewEx(MRI *mri, MRI *mri_mask, int orientation, int x,
                            int y, int z, int wsize, int fill_val, LTA *lta);
MRI *
MRIextractTalairachPlaneEx(MRI *mri_src, MRI *mri_dst, int orientation,
                           int x, int y, int z, int wsize, LTA *lta);
///////////////////////////////////////////////////////////////////////////////

#endif
