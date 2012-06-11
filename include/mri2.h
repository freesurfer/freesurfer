/**
 * @file  mri2.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: greve $
 *    $Date: 2012/06/11 17:42:47 $
 *    $Revision: 1.40 $
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


#ifndef MRI2_H
#define MRI2_H

#if defined(__cplusplus)
extern "C" {
#endif

#include "mri.h"
#include "mriTransform.h"
#include "mrisurf.h"
#include "connectcomp.h"

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
int MRIdimMismatch( const MRI *v1, const MRI *v2, int frameflag);
MATRIX *MRIcovarianceMatrix(MRI *mri, MRI *mask);
int MRIpca(MRI *D, MATRIX **pU, VECTOR **pS, MRI **pV, MRI *mask);
int WritePCAStats(char *fname, MATRIX *Spca);
int PrintPCAStats(FILE *fp, MATRIX *Spca);
MRI *MRIsqrt(MRI *invol, MRI *outvol);
double MRImaxAbsDiff(MRI *vol1, MRI *vol2,
                     int *cmax, int *rmax, int *smax, int *fmax);
MRI *MRImax(MRI *mri1, MRI *mri2, MRI *out);
MRI *MRImultiplyConst(MRI *src, double vconst, MRI *dst);
MRI *MRIaddConst(MRI *src, double vconst, MRI *dst);
int MRIvol2VolTkReg(MRI *mov, MRI *targ, MATRIX *Rtkreg,
		    int InterpCode, float param);
MRI *MRIvol2VolTLKernel(MRI *src, MRI *targ, MATRIX *Vt2s);
MRI *MRIvol2VolDelta(MRI *mov, MRI *targ, MATRIX *Rt2s);
MRI *MRIexp(MRI *mri, double a, double b, MRI *mask, MRI *out);
MRI *MRIsum(MRI *mri1, MRI *mri2, double a, double b, MRI *mask, MRI *out);
MRI *MRIvote(MRI *in, MRI *mask, MRI *vote);

#define VOX2VOXREGTYPE_FILE 0 /* Use specifed file */
#define VOX2VOXREGTYPE_FIND 1 /* Look for register.dat in movable MRI dir */
#define VOX2VOXREGTYPE_IDENTITY 2 /* Use MRItkRegMtx() */
int MRImakeVox2VoxReg(MRI* targ, MRI* mov,
                      int regtype, char* regname,
                      mriTransformRef* transform);
double MRIsum2All(MRI *mri);
MRI *MRIchecker(MRI *mri, MRI *checker);
MRI *MRIcrop(MRI *mri,int c1, int r1, int s1, int c2, int r2, int s2);
MRI *MRIuncrop(MRI *mri, MRI *crop, int c1, int r1, int s1, int c2, int r2, int s2);
MRI *MRIreverseSlices(MRI *in, MRI *out);
MRI *MRIcutEndSlices(MRI *mri, int ncut);
MRI *MRIsquare(MRI *in, MRI *mask, MRI *out);
MRI *MRIsquareRoot(MRI *in, MRI *mask, MRI *out);
int *MRIsegIdList(MRI *seg, int *nlist, int frame);
double *MRIsegDice(MRI *seg1, MRI *seg2, int *nsegs, int **segidlist);
MRI *MRIsegDiff(MRI *oldseg, MRI *newseg, int *DiffFlag);
MRI *MRIsegMergeDiff(MRI *oldseg, MRI *diff);
MRI *MRIhalfCosBias(MRI *in, double alpha, MRI *out);
int MRIvol2VolVSM(MRI *src, MRI *targ, MATRIX *Vt2s,
		  int InterpCode, float param, MRI *vsm);
int MRIvol2VolTkRegVSM(MRI *mov, MRI *targ, MATRIX *Rtkreg,
		       int InterpCode, float param, MRI *vsm);
MRI *MRIvol2surfVSM( const MRI *SrcVol,
                     const MATRIX *Rtk, const MRI_SURFACE *TrgSurf,
                     const MRI *vsm, int InterpMethod, MRI *SrcHitVol,
                     float ProjFrac, int ProjType, int nskip, 
		    MRI *TrgVol);
MRI *MRIsegBoundary(MRI *seg);
MRI *MRIsliceNo(MRI *in, MRI *out);
MRI *MRIindexNo(MRI *in, MRI *out);
MRI *MRIcrs(MRI *in, MRI *out);
int MRIsegFrameAvg(MRI *seg, int segid, MRI *mri, double *favg);
int MRIsegStats(MRI *seg, int segid, MRI *mri,int frame,
                float *min, float *max, float *range,
                float *mean, float *std);
int MRIsegStatsRobust(MRI *seg, int segid, MRI *mri,int frame,
		      float *min, float *max, float *range,
		      float *mean, float *std, float Pct);

MRI *MRImask_with_T2_and_aparc_aseg(MRI *mri_src, MRI *mri_dst, MRI *mri_T2, MRI *mri_aparc_aseg, float T2_thresh, int mm_from_exterior) ;
int *MRIsegmentationList(MRI *seg, int *pListLength);


#if defined(__cplusplus)
};
#endif

#endif
