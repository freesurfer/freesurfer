//
// talairachex.h
//
// Warning: Do not edit the following four lines.  CVS maintains them.
// Revision Author: $Author: fischl $
// Revision Date  : $Date: 2005/05/22 13:16:52 $
// Revision       : $Revision: 1.2 $
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
