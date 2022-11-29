/*
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


#ifdef X
#undef X
#endif

#ifndef MRI2_H
#define MRI2_H

#include "mri.h"
#include "mriTransform.h"
#include "mrisurf.h"
#include "connectcomp.h"

/*!
  \class class MRIvol2VolLTA
  Manages converting one volume to another via an LTA. All it needs is
  the mov and the LTA (it computes the output geometry from that in
  the LTA). Basic usage is to set the mov and set the LTA, then run
  vol2vol(); it will figure out the proper direction of the LTA.  If
  the LTA is not set but the target is, then uses a header
  registration. Handles a voxel-shift map (VSM). Note that the final
  result might be slightly different than mri_vol2vol because
  mri_vol2vol goes through several LTA type changes.
 */
class MRIvol2VolLTA{
public:
  LTA *lta=NULL;
  MRI *mov=NULL;
  MRI *targ=NULL;// only used for headerreg
  MRI *vsm=NULL;
  int InterpCode=SAMPLE_TRILINEAR;
  double sinchw=0;
  MRI *vol2vol(MRI *outvol);
  int ReadLTAorTarg(char *fname);
};


MRI *mri_load_bvolume(char *bfstem);
int  mri_save_as_bvolume(MRI *vol, char *stem, int svendian, int svtype);
MRI *mri_load_bvolume_frame(char *bfstem, int frameno);
int  mri_framepower(MRI *vol, float *framepower);
MRI *mri_binarize(MRI *vol, float thresh, const char *tail, int invert,
                  MRI *volbin, int *nover);
MRI *mri_rescale(MRI *invol, float min, float max, MRI *outvol);
int  mri_save_as_cor(MRI *vol,  char *cordir, int frame, int rescale);
int mri_minmax(MRI *vol, float *min, float *max);
MRI *mri_load_cor_as_float(char *cordir);
MRI *mri_load_wfile(char *wfile);
size_t mri_sizeof(MRI *vol);
MRI *mri_reshape(MRI *vol, int ncols, int nrows, int nslices, int nframes);
MRI *MRIreshape1d(MRI *src, MRI *trg);

int MRIfdr2vwth(MRI **vollist, int nvols, int *framelist, double fdr, int signid,
                int log10flag, MRI **masklist, double *vwth, MRI **ovollist);
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
int MRImostFreqNeighbor(MRI *mri, int c, int r, int s, int f, int delta);

#define VOX2VOXREGTYPE_FILE 0 /* Use specifed file */
#define VOX2VOXREGTYPE_FIND 1 /* Look for register.dat in movable MRI dir */
#define VOX2VOXREGTYPE_IDENTITY 2 /* Use MRItkRegMtx() */
int MRImakeVox2VoxReg(MRI* targ, MRI* mov,
                      int regtype, char* regname,
                      mriTransformRef* transform);
double MRIsum2All(MRI *mri);
MRI *MRIchecker(MRI *mri, MRI *checker);
MRI *MRIgrid(MRI *mri, int dc, int dr, int ds, float val, MRI *grid);
MRI *MRIcrop(MRI *mri,int c1, int r1, int s1, int c2, int r2, int s2);
MRI *MRIuncrop(MRI *mri, MRI *crop, int c1, int r1, int s1, int c2, int r2, int s2);
MRI *MRIreverseSlices(MRI *in, MRI *out);
MRI *MRIcutEndSlices(MRI *mri, int ncut);
MRI *MRIsquare(MRI *in, MRI *mask, MRI *out);
MRI *MRIsquareRoot(MRI *in, MRI *mask, MRI *out);
int *MRIsegIdList(MRI *seg, int *nlist, int frame);
int *MRIsegIdListNot0(MRI *seg, int *nsegs, int frame);
int *MRIsegIdListExclude0(MRI *seg, int *pnlist, int frame);
MRI *MRIbinarizeMatch(MRI *seg, int *MatchList, int nList, int frame, MRI *out);
double *MRIsegDice(MRI *seg1, MRI *seg2, int *nsegs, int **segidlist);
MRI *MRIsegDiff(MRI *oldseg, MRI *newseg, int *DiffFlag);
MRI *MRIsegMergeDiff(MRI *oldseg, MRI *diff);
MRI *MRIhalfCosBias(MRI *in, double alpha, MRI *out);
MRI *MRIvol2VolFill(MRI *src, MRI *mask, LTA *lta, int UpsampleFactor, int Average, MRI *outfill);
int MRIvol2VolVSM(MRI *src, MRI *targ, MATRIX *Vt2s,
		  int InterpCode, float param, MRI *vsm);
int MRIvol2VolTkRegVSM(MRI *mov, MRI *targ, MATRIX *Rtkreg,
		       int InterpCode, float param, MRI *vsm);
MRI *MRIvol2surfVSM( const MRI *SrcVol,
                     const MATRIX *Rtk, const MRI_SURFACE *TrgSurf,
                     const MRI *vsm, int InterpMethod, MRI *SrcHitVol,
                     float ProjFrac, int ProjType, int nskip, 
		    MRI *TrgVol);
MRI *MRImaskAndUpsample(MRI *src, MRI *mask, int UpsampleFactor, int nPad, int DoConserve, LTA **src2out);
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

MATRIX *BuildGTM0(MRI *seg, MRI *mask, double cFWHM, double rFWHM, double sFWHM, MATRIX *X);
MRI *MRIfisherTransform(MRI *rho, MRI *mask, MRI *out);
MRI *MRIhiresSeg(MRI *aseg, MRIS *lhw, MRIS *lhp, MRIS *rhw, MRIS *rhp, int USF, LTA **aseg2hrseg);
MRI *MRIpartialVolumeFraction(LTA *seg2vol, MRI *seg, double resmm, COLOR_TABLE *ct, MRI *pvf);
MRI *MRIpartialVolumeFractionAS(LTA *aseg2vol, MRI *aseg, MRIS *lhw, MRIS *lhp, 
				MRIS *rhw, MRIS *rhp, int USF, double resmm, COLOR_TABLE *ct, MRI *pvf);

int MRIcountMatches(const MRI *seg, const int MatchVal, const int frame, const MRI *mask);
MRI *MRIaddExtraCerebralCSF(MRI *seg, int nDil, MRI *out);
COLOR_TABLE *CTABpruneCTab(const COLOR_TABLE *ct0, MRI *seg);
MRI *MRIannot2CorticalSeg(MRI *seg, MRIS *lhw, MRIS *lhp, MRIS *rhw, MRIS *rhp, LTA *anat2seg, MRI *ctxseg);
MRI *MRIannot2CerebralWMSeg(MRI *seg, MRIS *lhw, MRIS *rhw, double DistThresh, LTA *anat2seg, MRI *wmseg);
MRI *MRIunsegmentCortex(MRI *seg, int lhmin, int lhmax, int rhmin, int rhmax, MRI *out);
MRI *MRIunsegmentWM(MRI *seg, MRIS *lhw, MRIS *rhw, int *segidlist, int nlist, LTA *anat2seg, MRI *wmseg);
MRI *MRIrelabelHypoHemi(MRI *seg, MRIS *lhw, MRIS *rhw, LTA *anat2seg, MRI *wmseg);
MRI *MRIrelabelNonWMHypos(MRI *seg0, int *segidlist, int nsegs, int *outsegidlist);
MRI *CTABcount2MRI(COLOR_TABLE *ct, MRI *seg);
MRI *MRIreorientLIA2RAS(MRI *mriA, MRI *mriB);

/*
  converts vertices to the index space
*/
void MRIConvertSurfaceVertexCoordinates(MRIS* mris, MRI* vol);

MATRIX *MRIvol2mat(MRI *vol, MRI *mask, int transposeFlag, MATRIX *M);
MRI *MRImat2vol(MATRIX *M, MRI *mask, int transposeFlag, MRI *vol);
MRI *MRImergeSegs(MRI *seg, int *seglist, int nsegs, int NewSegId, MRI *newseg);
MRI *MRImatchSegs(MRI *seg, int *seglist, int nsegs, int MaskId, MRI *mask);
HISTOGRAM *HISTOseg(MRI *seg, int segid, MRI *vol, double bmin, double bmax, double bdelta);
int QuadEulerCharChange(MRI *vol, MRI *mask, int c, int r, int s);
int QuadEulerCharChangeTest(int ForceFail);
int QuadEulerCharChangeCheckReorder(MRI *mri, const char *testname, int decExpected);

MRI *MRIfindBrightNonWM(MRI *mri_T1, MRI *mri_wm);
MRI *MRIzconcat(MRI *mri1, MRI *mri2, int nskip, MRI *out);
int MRIfillTriangle(MRI *vol, double p1[3], double p2[3], double p3[3], double dL, double FillVal);

/*!
  \fn class FixSubCortMassHA
  \brief This class "fixes" the SCM by removing stray voxels in
  hippocampus, amyg, and inferior lateral ventricle.  The SCM can be
  the wm.seg.mgz, wm.mgz, or filled.mgz.  The need for this stems from
  the fact that there is a lot of WM in the hippo (eg, alveus, fimria)
  and these often get segmented into the SCM where they create defects
  or make the surface distorted.  All the voxels in amyg and ILV are
  removed. Most of the voxels in hippo are removed except for those
  within nDilate of cortex, wm, and background. These are the voxels
  most likely to be in the wm of entorhinal cortex.
 */
class FixSubCortMassHA {
public:
  MRI *subcorticalmass=NULL; // wm.seg, wm, or filled
  MRI *aseg=NULL;  // aseg.presurf
  MRI *mask=NULL;  // a temporary mask
  int nDilate = 1;
  // These are the seg ids of structures outside the SCM, including cortex
  std::vector<int> OutsideSCMIds{0,165,258,259,24,3,42,17,18,53,54,85,5,44,31,63,16,7,8,46,47,15};
  int ExcludeOutsideSCM=0;
  MRI *scminsidemask=NULL;
  int nDilateInsideMask = 2;
  int FixSCM(void);
};

/*!
  \fn class BinarizeMRI
  \breif This class is used to help binarize MRIs. This can be done by
  simple thresholding or by matching values from a list. There are
  many bells and whistles, eg, inverting, applying to only certain
  frames, taking abs of the input before thesholding, constraining to
  a mask, allways zeroing edges.  A Voxel value must be
  thmin<=val<=thmax to be set; if invert=1, then thmin>val and
  val>thmax.
 */
class BinarizeMRI
{
public:
  int BinType=1; // 1=threshold, 2=match
  double thmin=0.5, thmax=INFINITY;
  std::vector<int> matchlist; // list of ints for when BinType=2
  double maskthmin = 0.5, maskthmax=INFINITY;
  int ZeroColEdges = 0;
  int ZeroRowEdges = 0;
  int ZeroSliceEdges = 0;
  double OnVal=1, OffVal=0;
  int fstart=-1, fend=-1;
  int invert=0;
  int nhits=0; // total number set to OnVal
  int mritype=-1; // default will be MRI_FLOAT
  int DoAbs=0;
  int m_debug=0;
  FILE *m_debugfp=stdout;
  int qualifies(const MRI *invol, const MRI *mask,  int c,  int r,  int s,  int f);
  MRI *binarize(MRI *invol, const MRI *mask, MRI *outvol);
  int dump(MRI *invol, const MRI *mask, MRI *outvol);
  // add: fdr, match-from-{ctab,ctx,gm,wm,subctx}, bincol, frame ops (sum,min,max,and,or)
};

/*!
  \fn class DEMorphBinVol
  \brief Performs dilation (morphtype=1) and erosion (morphtype=2)
  type morphological operations on binary MRI volumes. Combining topo
  with min number of neighbors can allow the extent of
  dilation/erosion to be controlled better. See class BinarizeMRI 
  for many options for binarizing before calling these methods.
 */
class DEMorphBinVol
{
public:
  int morphtype=1; //1=dil, 2=erode
  int topo=1; // 1=face neighors, 2=face+edges neighbors, 3=face+edge+corner neighbors
  int nnbrsthresh=0; // min number of tagged neighbors must exceed this number
  int nmorph=0; //number of iterations, 0=do nothing
  int frame=0; // use this frame from the input
  int nchangestot=0; // set to the number of voxels changed during morph
  int mm_mritype=-1; // -1 = make output same type as input
  int mm_debug=1;
  FILE *mm_debugfp=stdout;
  int check(MRI *invol, MRI *mask, MRI *outvol);
  int dump(MRI *binvol);
  MRI *copyFrame(MRI *invol, int frameno);
  MRI *morph(MRI *binvol); // basic morph of a binary volumes
  // These two will binarize for you, then run the morph
  MRI *morph(MRI *invol, std::vector<int> matchlist, MRI *mask, MRI *outvol);
  MRI *morph(MRI *invol, double thmin, double thmax, MRI *mask, MRI *outvol);
};

/*!
  \fn class SCMstopMask
  \brief This class supplies methods to create a "stop mask" for white
  surface placement. The stopmask is used when CBV is searching for
  peaks in the gradient. By stopping the search, problematic voxels
  are avoided (eg, the occip horn of the LVs) and the surface does not
  wander out into silly areas. The mask can be created based on
  several criteria: (1) Lateral ventricles/choroid plexus (needs
  aseg.presurf), (2) WMSAs (77,78,79, needs aseg.presurf); WMSAs can
  be tricky because cortex is sometimes mislabled as WMSA; therefore,
  the ability to exclude WMSA voxels near cortex has been included (it
  is a bit of a hack). (3) An active edit to the filled (needs
  filled.auto.mgz and filled.mgz). (4) an active edit to the wm.mgz
  (5) an active edit to the brain.finalsurfs. Ideally, one would only
  be editing the filled.mgz; the others are there for backwards
  compatibility and because the long does not currently propagate the
  filled. The output is actually a segmentation with colortable; the
  segid indicates which criterion was used to set the voxel.
 */
class SCMstopMask {
public:
  MRI *bfs=NULL; // brain.finalsurfs.mgz
  MRI *aseg=NULL; // aseg.presurf.mgz
  MRI *wm=NULL; // wm.mgz
  MRI *filledauto=NULL, *filled=NULL; // filled.auto.mgz and filled.mgz
  int DoBFS255 = 1; // add voxels in BFS=255 (these have been manually edited)
  int DoFilled=1; // add filled voxels that are > 0 and have been edited
  int DoLV=1; // add lateral ventricle and choroid plexus voxels (add ILV too?)
  int DoWM255=1; // add edited WM voxels=255
  int DoWMSA=1; // add wm hypointensites (77,78,79)
  double WMSAErodeMM=1; // erode any WMSAs that are adjacent to cortex by this amount in mm
  int nWMSAErode; // number of voxels to erode WMSAs, computed later
  // Add putamen?
  // Could find places in aseg that have been edited to be part of the SCM and then force
  // them into the stopmask
  MRI *getmask(void);
  int nhitslist[6]; // to keep track of number of voxels hit
};


#endif
