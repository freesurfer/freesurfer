/**
 * @file  mri2.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 02:09:00 $
 *    $Revision: 1.16 $
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


#ifndef MRI2_H
#define MRI2_H

#include "mri.h"
#include "mriTransform.h"

MRI *mri_load_bvolume(char *bfstem);
int  mri_save_as_bvolume(MRI *vol, char *stem, int svendian, int svtype);
MRI *mri_load_bvolume_frame(char *bfstem, int frameno);
int  mri_framepower(MRI *vol, float *framepower);
MRI *mri_binarize(MRI *vol, float thresh, char *tail, int invert,
                  MRI *volbin, int *nover);
MRI *mri_rescale(MRI *invol, float min, float max, MRI *outvol);
int  mri_save_as_cor(MRI *vol,  char *cordir, int frame, int rescale);
int mri_minmax(MRI *vol, float *min, float *max);
MRI *mri_load_cor_as_float(char *cordir);
MRI *mri_load_wfile(char *wfile);
size_t mri_sizeof(MRI *vol);
MRI *mri_reshape(MRI *vol, int ncols, int nrows, int nslices, int nframes);
int MRIfdr2vwth(MRI *vol, int frame, double fdr, int signid,
                int log10flag, MRI *mask, double *vwth, MRI *ovol);
int MRIdimMismatch(MRI *v1, MRI *v2, int frameflag);
MATRIX *MRIcovarianceMatrix(MRI *mri, MRI *mask);
int MRIpca(MRI *D, MATRIX **pU, VECTOR **pS, MRI **pV, MRI *mask);
int WritePCAStats(char *fname, MATRIX *Spca);
int PrintPCAStats(FILE *fp, MATRIX *Spca);
MRI *MRIsqrt(MRI *invol, MRI *outvol);
double MRImaxAbsDiff(MRI *vol1, MRI *vol2,
                     int *cmax, int *rmax, int *smax, int *fmax);
MRI *MRImax(MRI *mri1, MRI *mri2, MRI *out);
MRI *MRImultiplyConst(MRI *src, double vconst, MRI *dst);
MRI *MRIvol2VolTLKernel(MRI *src, MRI *targ, MATRIX *Vt2s);
MRI *MRIexp(MRI *mri, double a, double b, MRI *mask, MRI *out);
MRI *MRIsum(MRI *mri1, MRI *mri2, double a, double b, MRI *mask, MRI *out);

#define VOX2VOXREGTYPE_FILE 0 /* Use specifed file */
#define VOX2VOXREGTYPE_FIND 1 /* Look for register.dat in movable MRI dir */
#define VOX2VOXREGTYPE_IDENTITY 2 /* Use MRItkRegMtx() */
int MRImakeVox2VoxReg(MRI* targ, MRI* mov,
                      int regtype, char* regname,
                      mriTransformRef* transform);
#endif
