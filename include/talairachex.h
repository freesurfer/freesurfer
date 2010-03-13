/**
 * @file  talairachex.h
 * @brief new talairach related routines with Ex
 *
 * takes lta as the talairach transform (use LTAreadEx routine)
 * doesn't rely on COR volume type
 */
/*
 * Original Author: Y. Tosa
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2010/03/13 01:32:40 $
 *    $Revision: 1.5 $
 *
 * Copyright (C) 2002-2010,
 * The General Hospital Corporation (Boston, MA). 
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 *
 */

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
void TransformWithMatrix( const MATRIX *mat,
			  const double x,
			  const double y,
			  const double z,
			  double *px,
			  double *py,
			  double *pz );

int
ModifyTalairachCRAS(MRI *mri_tal, const LTA *lta);
int
MRIvoxelToTalairachEx(MRI *mri_src, double xv, double yv, double zv,
                      double *pxt, double *pyt, double *pzt, const LTA *lta);
int
MRItalairachToVoxelEx(MRI *mri_dst, double xt, double yt, double zt,
                      double *pxv, double *pyv, double *pzv, const LTA *lta);
int
MRItalairachVoxelToWorldEx(MRI *mri_dst,
                           double xtv, double ytv, double ztv,
                           double *pxw, double *pyw, double *pzw,
                           const LTA *lta);
int
MRIvoxelToTalairachVoxelEx(MRI *mri_src,
                           double xv, double yv, double zv,
                           double *pxt, double *pyt, double *pzt,
                           const LTA *lta);
int
MRItalairachVoxelToVoxelEx(MRI *mri_dst,
                           double xv, double yv, double zv,
                           double *pxnv, double *pynv, double *pznv,
                           const LTA *lta) ;
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
