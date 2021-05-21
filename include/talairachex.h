/**
 * @brief new talairach related routines with Ex
 *
 * takes lta as the talairach transform (use LTAreadEx routine)
 * doesn't rely on COR volume type
 */
/*
 * Original Author: Y. Tosa
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
